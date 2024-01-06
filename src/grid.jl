"""
function GridInfo(N_real::NTuple{3, Int}, w::NTuple{3, Int}, periodicity::NTuple{3, Bool}, pad::NTuple{3, Int}, L::NTuple{3, T}) where{T<:Union{Float32, Float64}}

    This function creates a GridBox object with two mutable Array.
    
    N_real: the number of real grid points in each direction
    w: the number of grid points in each direction for the window function cutoff
    periodicity: a tuple of three boolean values indicating whether the grid is periodic in each direction, to decide padding or not
    extra_pad: the number of extra grid points in each direction for padding
    L: the length of the box in each direction
    h: the grid spacing in each direction, the first point is at [h/2, 3h/2, ..., L - h/2]
    N_image: the number of grid points in each direction for the image grid, given by tuple([2 * w[i] + 2 for i in 1:3]...) .+ N_real
    N_pad: the numer of padded grid points in each direction, given by tuple([2 * pad[i] * (periodicity[i] == false) for i in 1:3]...) .+ N_image .+ N_real
"""
function GridInfo(N_real::NTuple{3, Int}, w::NTuple{3, Int}, periodicity::NTuple{3, Bool}, extra_pad::NTuple{3, Int}, L::NTuple{3, T}) where{T<:Union{Float32, Float64}}
    N_image = tuple([2 * w[i] + 2 for i in 1:3]...) .+ N_real
    image = tuple([w[i] + 1 for i in 1:3]...)
    pad = tuple([(extra_pad[i] + w[i] + 1) * (periodicity[i] == false) for i in 1:3]...)
    N_pad = tuple([2 * pad[i] for i in 1:3]...) .+ N_real
    h = tuple([L[i] / (N_real[i]) for i in 1:3]...)

    trans_info = TransInfo(N_real, periodicity, image, pad)

    return GridInfo{T}(N_real, w, periodicity, image, pad, L, h, N_image, N_pad, trans_info)
end

function GridBox(grid_info::GridInfo{T}) where{T<:Union{Float32, Float64}}
    pad_grid = zeros(T, grid_info.N_pad...)
    image_grid = zeros(T, grid_info.N_image...)
    return GridBox{T}(pad_grid, image_grid)
end

"""
function TransInfo(N_real::NTuple{3, Int}, periodicity::NTuple{3, Bool}, image::NTuple{3, Int}, pad::NTuple{3, Int})

    this function is used the generate the transformation information for the image grid and the pad grid
    the space will be divided into 27 regions, each region has a transformation information in form of Vector of Tuples given by ((ix, iy, iz), (px, py, pz), (sx, sy, sz))

    When doing grid transformation, use it as:
    for i in 1:sx
        for j in 1:sy
            for k in 1:sz
                pad_grid[px + i, py + j, pz + k] += image_grid[ix + i, iy + j, iz + k]
            end
        end
    end
"""
function TransInfo(N_real::NTuple{3, Int}, periodicity::NTuple{3, Bool}, image::NTuple{3, Int}, pad::NTuple{3, Int})

    # the first tuple for image grid, the second tuple for pad grid, the third tuple for size
    trans_info = Vector{Tuple{NTuple{3, Int}, NTuple{3, Int}, NTuple{3, Int}}}()

    relative_pos = (-1, 0, 1)

    start_image = [(0, image[i], image[i] + N_real[i]) for i in 1:3]
    start_pad = [(pad[i] - image[i], pad[i], pad[i] + N_real[i]) for i in 1:3]
    size_trans = [(image[i], N_real[i], image[i]) for i in 1:3]
    
    for cart in CartesianIndices((3, 3, 3))
        i, j, k = cart[1], cart[2], cart[3]
        ix, iy, iz = start_image[1][i], start_image[2][j], start_image[3][k]

        px = start_pad[1][i] - relative_pos[i] * N_real[1] * (periodicity[1] == true)
        py = start_pad[2][j] - relative_pos[j] * N_real[2] * (periodicity[2] == true)
        pz = start_pad[3][k] - relative_pos[k] * N_real[3] * (periodicity[3] == true)

        sx, sy, sz = size_trans[1][i], size_trans[2][j], size_trans[3][k]

        push!(trans_info, ((ix, iy, iz), (px, py, pz), (sx, sy, sz)))
    end

    return trans_info
end

"""
function grid_revise_image!(gridbox::GridBox{T}) where{T}

    set all values of gridbox.image_grid to 0
"""
function grid_revise_image!(gridbox::GridBox{T}) where{T}

    for i in eachindex(gridbox.image_grid)
        gridbox.image_grid[i] = zero(T)
    end

    return nothing
end

"""
function grid_revise_pad!(gridbox::GridBox{T}) where{T}

    set all values of gridbox.pad_grid to 0
"""
function grid_revise_pad!(gridbox::GridBox{T}) where{T}

    for i in eachindex(gridbox.pad_grid)
        gridbox.pad_grid[i] = zero(T)
    end

    return nothing
end

@inbounds function grid_image2pad!(gridbox::GridBox{T}, gridinfo::GridInfo{T}) where T

    grid_revise_pad!(gridbox)
    trans_info = gridinfo.trans_info
    for (iid, pid, sid) in trans_info
        (ix, iy, iz), (px, py, pz), (sx, sy, sz) = iid, pid, sid
        @turbo for i in 1:sx
            for j in 1:sy
                for k in 1:sz
                    gridbox.pad_grid[px + i, py + j, pz + k] += gridbox.image_grid[ix + i, iy + j, iz + k]
                end
            end
        end
    end

    return gridbox
end

@inbounds function grid_pad2image!(gridbox::GridBox{T}, gridinfo::GridInfo{T}) where T

    grid_revise_image!(gridbox)
    trans_info = gridinfo.trans_info
    for (iid, pid, sid) in trans_info
        (ix, iy, iz), (px, py, pz), (sx, sy, sz) = iid, pid, sid
        @turbo for i in 1:sx
            for j in 1:sy
                for k in 1:sz
                    gridbox.image_grid[ix + i, iy + j, iz + k] = gridbox.pad_grid[px + i, py + j, pz + k]
                end
            end
        end
    end

    return gridbox
end