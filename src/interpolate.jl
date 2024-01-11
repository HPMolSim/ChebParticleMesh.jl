function interpolate_single!(q::T, pos::NTuple{N, T}, gridinfo::GridInfo{N, T}, gridbox::GridBox{N, T}, chebcoefs::NTuple{N, ChebCoef{T}}) where{N, T}

    cheb_value = gridbox.cheb_value
    near_id_image = image_grid_id(pos, gridinfo)
    near_pos_image = image_grid_pos(near_id_image, gridinfo)
    for i in 1:N
        dx = pos[i] - near_pos_image[i]
        pwcheb_eval!(dx, cheb_value[i], chebcoefs[i])
    end

    for id in Iterators.product([- gridinfo.w[i] : gridinfo.w[i] for i in 1:N]...)
        gridbox.image_grid[(near_id_image.id .+ id)...] += q * prod(cheb_value[i][id[i] + gridinfo.w[i] + 1] for i in 1:N)
    end

    return nothing
end

function interpolate!(qs::Vector{T}, poses::Vector{NTuple{N, T}}, gridinfo::GridInfo{N, T}, gridbox::GridBox{N, T}, chebcoefs::NTuple{N, ChebCoef{T}}) where{N, T}

    @assert length(qs) == length(poses)
    grid_revise_image!(gridbox)

    for i in 1:length(qs)
        interpolate_single!(qs[i], poses[i], gridinfo, gridbox, chebcoefs)
    end

    return nothing
end

function interpolate_single_direct!(q::T, pos::NTuple{N, T}, gridinfo::GridInfo{N, T}, gridbox::GridBox{N, T}, chebcoefs::NTuple{N, ChebCoef{T}}) where{N, T}

    cheb_value = gridbox.cheb_value
    near_id_image = image_grid_id(pos, gridinfo)
    near_pos_image = image_grid_pos(near_id_image, gridinfo)
    for i in 1:N
        dx = pos[i] - near_pos_image[i]
        f_eval!(dx, cheb_value[i], chebcoefs[i])
    end

    for id in Iterators.product([- gridinfo.w[i] : gridinfo.w[i] for i in 1:N]...)
        gridbox.image_grid[(near_id_image.id .+ id)...] += q * prod(cheb_value[i][id[i] + gridinfo.w[i] + 1] for i in 1:N)
    end

    return nothing
end

function interpolate_direct!(qs::Vector{T}, poses::Vector{NTuple{N, T}}, gridinfo::GridInfo{N, T}, gridbox::GridBox{N, T}, chebcoefs::NTuple{N, ChebCoef{T}}) where{N, T}

    @assert length(qs) == length(poses)
    grid_revise_image!(gridbox)

    for i in 1:length(qs)
        interpolate_single_direct!(qs[i], poses[i], gridinfo, gridbox, chebcoefs)
    end

    return nothing
end