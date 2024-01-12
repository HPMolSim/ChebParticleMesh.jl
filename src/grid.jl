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

    index_list = [mod1.(1 + pad[i] - image[i]:pad[i] + N_real[i] + image[i], N_pad[i]) for i in 1:N]
    iter_list = Iterators.product([- w[i] : w[i] for i in 1:N]...)

    return GridInfo{N, T}(N_real, w, periodicity, image, pad, L, h, N_image, N_pad, index_list, iter_list, k)
end

function GridBox(grid_info::GridInfo{N, T}) where{N, T<:Union{Float32, Float64}}
    pad_grid = zeros(Complex{T}, grid_info.N_pad...)
    
    image_grid = view(pad_grid, grid_info.index_list...)

    cheb_value = [zeros(T, 2 * grid_info.w[i] + 1) for i in 1:N]
    return GridBox{N, T}(pad_grid, image_grid, cheb_value)
end


"""
function grid_revise_pad!(gridbox::GridBox{N, T}) where{N, T}

    set all values of gridbox.pad_grid to 0
"""
function grid_revise_pad!(gridbox::GridBox{N, T}) where{N, T}

    for i in eachindex(gridbox.pad_grid)
        gridbox.pad_grid[i] = zero(Complex{T})
    end

    return nothing
end