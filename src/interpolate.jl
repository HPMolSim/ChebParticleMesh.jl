function interpolate_single!(q::T, pos::NTuple{N, T}, gridinfo::GridInfo{N, T}, gridbox::GridBox{N, T}, chebcoefs::NTuple{N, ChebCoef{T}}) where{N, T}

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
        value = Complex{T}(q * prod(cheb_value[1][i[j] + gridinfo.w[1] + 1] for j in 1:N))
        revise_pad!(gridbox.pad_grid, idl, image_id, value)
    end

    return nothing
end

function interpolate!(qs::Vector{T}, poses::Vector{NTuple{N, T}}, gridinfo::GridInfo{N, T}, gridbox::GridBox{N, T}, chebcoefs::NTuple{N, ChebCoef{T}}) where{N, T}

    @assert length(qs) == length(poses)
    grid_revise_pad!(gridbox)

    for i in 1:length(qs)
        interpolate_single!(qs[i], poses[i], gridinfo, gridbox, chebcoefs)
    end

    return nothing
end

function interpolate_single_direct!(q::T, pos::NTuple{N, T}, gridinfo::GridInfo{N, T}, gridbox::GridBox{N, T}, chebcoefs::NTuple{N, ChebCoef{T}}) where{N, T}

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
        value = Complex{T}(q * prod(cheb_value[1][i[j] + gridinfo.w[1] + 1] for j in 1:N))
        revise_pad!(gridbox.pad_grid, idl, image_id, value)
    end

    return nothing
end

function interpolate_direct!(qs::Vector{T}, poses::Vector{NTuple{N, T}}, gridinfo::GridInfo{N, T}, gridbox::GridBox{N, T}, chebcoefs::NTuple{N, ChebCoef{T}}) where{N, T}

    @assert length(qs) == length(poses)
    grid_revise_pad!(gridbox)

    for i in 1:length(qs)
        interpolate_single_direct!(qs[i], poses[i], gridinfo, gridbox, chebcoefs)
    end

    return nothing
end