"""
function GridInfo(N_real::NTuple{N, Int}, w::NTuple{N, Int}, periodicity::NTuple{N, Bool}, extra_pad::NTuple{N, Int}, L::NTuple{N, T}) where{N, T<:Union{Float32, Float64}}

    This function creates a GridBox object with two mutable Array.
    
    N_real: the number of real grid points in each direction
    w: the number of grid points in each direction for the window function cutoff
    periodicity: a tuple of three boolean values indicating whether the grid is periodic in each direction, to decide padding or not
    extra_pad: the number of extra grid points in each direction for padding
    L: the length of the box in each direction
    h: the grid spacing in each direction, the first point is at [h/2, 3h/2, ..., L - h/2]
    N_image: the number of grid points in each direction for the image grid, given by tuple([2 * w[i] + 2 for i in 1:N]...) .+ N_real
    N_pad: the numer of padded grid points in each direction, given by tuple([2 * pad[i] * (periodicity[i] == false) for i in 1:N]...) .+ N_image .+ N_real
"""
function GridInfo(N_real::NTuple{N, Int}, w::NTuple{N, Int}, periodicity::NTuple{N, Bool}, extra_pad::NTuple{N, Int}, L::NTuple{N, T}) where{N, T<:Union{Float32, Float64}}
    N_image = tuple([2 * w[i] + 2 for i in 1:N]...) .+ N_real
    image = tuple([w[i] + 1 for i in 1:N]...)
    pad = tuple([(extra_pad[i] + w[i] + 1) * (periodicity[i] == false) for i in 1:N]...)
    N_pad = tuple([2 * pad[i] for i in 1:N]...) .+ N_real
    h = tuple([L[i] / (N_real[i]) for i in 1:N]...)

    trans_info = TransInfo(N_real, periodicity, image, pad)

    k = Vector{Array{T, 1}}()

    for i in 1:N
        ki = [2π * j / (L[i] + h[i] * 2 * pad[i]) for j in 0:(N_pad[i] - 1)]

        for j in 1:N_pad[i]
            if j - 1 > ceil(N_pad[i] / 2)
                ki[j] -= 2π * N_pad[i] / (L[i] + h[i] * 2 * pad[i])
            end
        end

        push!(k, ki)
    end

    return GridInfo{N, T}(N_real, w, periodicity, image, pad, L, h, N_image, N_pad, trans_info, k)
end

function GridBox(grid_info::GridInfo{N, T}) where{N, T<:Union{Float32, Float64}}
    pad_grid = zeros(T, grid_info.N_pad...)
    image_grid = zeros(T, grid_info.N_image...)
    cheb_value = [zeros(T, 2 * grid_info.w[i] + 1) for i in 1:N]
    return GridBox{N, T}(pad_grid, image_grid, cheb_value)
end

"""
function TransInfo(N_real::NTuple{N, Int}, periodicity::NTuple{N, Bool}, image::NTuple{N, Int}, pad::NTuple{N, Int}) where{N}

    this function is used the generate the transformation information for the image grid and the pad grid
    the space will be divided into 27 regions, each region has a transformation information in form of Vector of Tuples given by (image_id, pad_id, size)

    When doing grid transformation, use it as:
    for i in 1:sx
        for j in 1:sy
            for k in 1:sz
                pad_grid[px + i, py + j, pz + k] += image_grid[ix + i, iy + j, iz + k]
            end
        end
    end
"""
function TransInfo(N_real::NTuple{N, Int}, periodicity::NTuple{N, Bool}, image::NTuple{N, Int}, pad::NTuple{N, Int}) where{N}

    # the first tuple for image grid, the second tuple for pad grid, the third tuple for size
    trans_info = Vector{Tuple{NTuple{N, UnitRange{Int}}, NTuple{N, UnitRange{Int}}}}()

    relative_pos = (-1, 0, 1)

    start_image = [(0, image[i], image[i] + N_real[i]) for i in 1:N]
    start_pad = [(pad[i] - image[i], pad[i], pad[i] + N_real[i]) for i in 1:N]
    size_trans = [(image[i], N_real[i], image[i]) for i in 1:N]
    
    for cart in CartesianIndices(tuple([3 for i in 1:N]...))

        iid = [start_image[i][cart[i]] for i in 1:N]
        pid = [start_pad[i][cart[i]] - relative_pos[cart[i]] * N_real[i] * (periodicity[i] == true) for i in 1:N]
        sizeinfo = [size_trans[i][cart[i]] for i in 1:N]

        image_region = tuple([iid[i] + 1 : iid[i] + sizeinfo[i] for i in 1:N]...)
        pad_region = tuple([pid[i] + 1 : pid[i] + sizeinfo[i] for i in 1:N]...)

        push!(trans_info, (image_region, pad_region))
    end

    return trans_info
end

"""
function grid_revise_image!(gridbox::GridBox{N, T}) where{N, T}

    set all values of gridbox.image_grid to 0
"""
function grid_revise_image!(gridbox::GridBox{N, T}) where{N, T}

    for i in eachindex(gridbox.image_grid)
        gridbox.image_grid[i] = zero(T)
    end

    return nothing
end

"""
function grid_revise_pad!(gridbox::GridBox{N, T}) where{N, T}

    set all values of gridbox.pad_grid to 0
"""
function grid_revise_pad!(gridbox::GridBox{N, T}) where{N, T}

    for i in eachindex(gridbox.pad_grid)
        gridbox.pad_grid[i] = zero(T)
    end

    return nothing
end

@inbounds function grid_image2pad!(gridbox::GridBox{2, T}, gridinfo::GridInfo{2, T}) where {T}

    grid_revise_pad!(gridbox)
    trans_info = gridinfo.trans_info
    for (image_region, pad_region) in trans_info
        @turbo for i in 1:length(image_region[1])
            for j in 1:length(image_region[2])
                gridbox.pad_grid[pad_region[1][i], pad_region[2][j]] += gridbox.image_grid[image_region[1][i], image_region[2][j]]
            end
        end
    end

    return gridbox
end

@inbounds function grid_pad2image!(gridbox::GridBox{2, T}, gridinfo::GridInfo{2, T}) where {T}

    grid_revise_pad!(gridbox)
    trans_info = gridinfo.trans_info
    for (image_region, pad_region) in trans_info
        @turbo for i in 1:length(image_region[1])
            for j in 1:length(image_region[2])
                gridbox.image_grid[image_region[1][i], image_region[2][j]] += gridbox.pad_grid[pad_region[1][i], pad_region[2][j]]
            end
        end
    end

    return gridbox
end

@inbounds function grid_image2pad!(gridbox::GridBox{3, T}, gridinfo::GridInfo{3, T}) where {T}

    grid_revise_pad!(gridbox)
    trans_info = gridinfo.trans_info
    for (image_region, pad_region) in trans_info
        @turbo for i in 1:length(image_region[1])
            for j in 1:length(image_region[2])
                for k in 1:length(image_region[3])
                    gridbox.pad_grid[pad_region[1][i], pad_region[2][j], pad_region[3][k]] += gridbox.image_grid[image_region[1][i], image_region[2][j], image_region[3][k]]
                end
            end
        end
    end

    return gridbox
end

@inbounds function grid_pad2image!(gridbox::GridBox{3, T}, gridinfo::GridInfo{3, T}) where {T}

    grid_revise_pad!(gridbox)
    trans_info = gridinfo.trans_info
    for (image_region, pad_region) in trans_info
        @turbo for i in 1:length(image_region[1])
            for j in 1:length(image_region[2])
                for k in 1:length(image_region[3])
                    gridbox.image_grid[image_region[1][i], image_region[2][j], image_region[3][k]] += gridbox.pad_grid[pad_region[1][i], pad_region[2][j], pad_region[3][k]]
                end
            end
        end
    end

    return gridbox
end

@inbounds function grid_image2pad!(gridbox::GridBox{N, T}, gridinfo::GridInfo{N, T}) where {N, T}

    grid_revise_pad!(gridbox)
    trans_info = gridinfo.trans_info
    for (image_region, pad_region) in trans_info
        (@view gridbox.pad_grid[pad_region...]) .+= (@view gridbox.image_grid[image_region...])
    end

    return gridbox
end

@inbounds function grid_pad2image!(gridbox::GridBox{N, T}, gridinfo::GridInfo{N, T}) where {N, T}

    grid_revise_image!(gridbox)
    trans_info = gridinfo.trans_info
    for (image_region, pad_region) in trans_info
        (@view gridbox.image_grid[image_region...]) .+= (@view gridbox.pad_grid[pad_region...])
    end

    return gridbox
end