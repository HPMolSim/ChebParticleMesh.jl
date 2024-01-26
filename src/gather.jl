@inbounds function gather_single(q::T, pos::NTuple{N, T}, gridinfo::GridInfo{N, T}, gridbox::GridBox{N, T}, chebcoefs::NTuple{N, ChebCoef{T}}) where{N, T}
    
    potential_i = zero(T)

    cheb_value = gridbox.cheb_value
    idl = gridinfo.index_list

    near_id_image = image_grid_id(pos, gridinfo)
    near_pos_image = image_grid_pos(near_id_image, gridinfo)
    for i in 1:N
        dx = pos[i] - near_pos_image[i]
        pwcheb_eval!(dx, cheb_value[i], chebcoefs[i])
    end

    for i in gridinfo.iter_list
        image_id = near_id_image.id .+ i
        potential_i += real(gridbox.pad_grid[ntuple(j -> idl[j][image_id[j]], N)...]) * prod(cheb_value[j][i[j] + gridinfo.w[j] + 1] for j in 1:N)
    end

    return q * 4π * prod(gridinfo.h) * potential_i
end

@inbounds function gather(qs::Vector{T}, poses::Vector{NTuple{N, T}}, gridinfo::GridInfo{N, T}, gridbox::GridBox{N, T}, chebcoefs::NTuple{N, ChebCoef{T}}) where{N, T}

    @assert length(qs) == length(poses)

    potential = zero(T)
    for i in 1:length(qs)
        potential += gather_single(qs[i], poses[i], gridinfo, gridbox, chebcoefs)
    end

    return potential
end

@inbounds function gather_single_direct(q::T, pos::NTuple{N, T}, gridinfo::GridInfo{N, T}, gridbox::GridBox{N, T}, chebcoefs::NTuple{N, ChebCoef{T}}) where{N, T}
    
    potential_i = zero(T)

    cheb_value = gridbox.cheb_value
    idl = gridinfo.index_list

    near_id_image = image_grid_id(pos, gridinfo)
    near_pos_image = image_grid_pos(near_id_image, gridinfo)
    for i in 1:N
        dx = pos[i] - near_pos_image[i]
        f_eval!(dx, cheb_value[i], chebcoefs[i])
    end

    for i in gridinfo.iter_list
        image_id = near_id_image.id .+ i
        potential_i += real(gridbox.pad_grid[ntuple(j -> idl[j][image_id[j]], N)...]) * prod(cheb_value[j][i[j] + gridinfo.w[j] + 1] for j in 1:N)
    end

    return q * 4π * prod(gridinfo.h) * potential_i
end

@inbounds function gather_direct(qs::Vector{T}, poses::Vector{NTuple{N, T}}, gridinfo::GridInfo{N, T}, gridbox::GridBox{N, T}, chebcoefs::NTuple{N, ChebCoef{T}}) where{N, T}

    @assert length(qs) == length(poses)

    potential = zero(T)
    for i in 1:length(qs)
        potential += gather_single_direct(qs[i], poses[i], gridinfo, gridbox, chebcoefs)
    end

    return potential
end

function gather_single_naive(q::T, pos::NTuple{N, T}, gridinfo::GridInfo{N, T}, gridbox::GridBox{N, T}, chebcoefs::NTuple{N, ChebCoef{T}}) where{N, T}
    
    potential_i = zero(T)

    near_id_image = image_grid_id(pos, gridinfo)

    for d_id in Iterators.product([(-gridinfo.image[d]):(gridinfo.image[d]) for d in 1:N]...)
        image_id = near_id_image.id .+ d_id
        pos_image = image_grid_pos(ImageIndex(image_id), gridinfo)
        d_pos = pos .- pos_image
        potential_i += real(gridbox.image_grid[image_id...]) * prod(chebcoefs[d].f(d_pos[d]) for d in 1:N)
    end

    return q * 4π * prod(gridinfo.h) * potential_i
end