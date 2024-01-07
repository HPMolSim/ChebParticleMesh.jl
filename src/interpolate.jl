function interpolate_single!(q::T, pos::Point{3, T}, gridinfo::GridInfo{T}, gridbox::GridBox{T}, chebcoefs::NTuple{3, ChebCoef{T}}) where{T}

    cheb_value = gridbox.cheb_value
    near_id_image = image_grid_id(pos, gridinfo)
    near_pos_image = image_grid_pos(near_id_image, gridinfo)
    for i in 1:3
        dx = pos[i] - near_pos_image[i]
        pwcheb_eval!(dx, cheb_value[i], chebcoefs[i])
    end

    @turbo for i in - gridinfo.w[1] : gridinfo.w[1]
        for j in - gridinfo.w[2] : gridinfo.w[2]
            for k in - gridinfo.w[3] : gridinfo.w[3]
                gridbox.image_grid[
                    near_id_image.idx + i, 
                    near_id_image.idy + j, 
                    near_id_image.idz + k] += 
                    q * 
                    cheb_value[1][i + gridinfo.w[1] + 1] * 
                    cheb_value[2][j + gridinfo.w[2] + 1] * 
                    cheb_value[3][k + gridinfo.w[3] + 1]
            end
        end
    end

    return nothing
end

function interpolate!(qs::Vector{T}, poses::Vector{Point{3, T}}, gridinfo::GridInfo{T}, gridbox::GridBox{T}, chebcoefs::NTuple{3, ChebCoef{T}}) where{T}

    @assert length(qs) == length(poses)
    grid_revise_image!(gridbox)

    for i in 1:length(qs)
        interpolate_single!(qs[i], poses[i], gridinfo, gridbox, chebcoefs)
    end

    return nothing
end

function interpolate_single_direct!(q::T, pos::Point{3, T}, gridinfo::GridInfo{T}, gridbox::GridBox{T}, chebcoefs::NTuple{3, ChebCoef{T}}) where{T}

    cheb_value = gridbox.cheb_value
    near_id_image = image_grid_id(pos, gridinfo)
    near_pos_image = image_grid_pos(near_id_image, gridinfo)
    for i in 1:3
        dx = pos[i] - near_pos_image[i]
        f_eval!(dx, cheb_value[i], chebcoefs[i])
    end

    @turbo for i in - gridinfo.w[1] : gridinfo.w[1]
        for j in - gridinfo.w[2] : gridinfo.w[2]
            for k in - gridinfo.w[3] : gridinfo.w[3]
                gridbox.image_grid[
                    near_id_image.idx + i, 
                    near_id_image.idy + j, 
                    near_id_image.idz + k] += 
                    q * 
                    cheb_value[1][i + gridinfo.w[1] + 1] * 
                    cheb_value[2][j + gridinfo.w[2] + 1] * 
                    cheb_value[3][k + gridinfo.w[3] + 1]
            end
        end
    end

    return nothing
end

function interpolate_direct!(qs::Vector{T}, poses::Vector{Point{3, T}}, gridinfo::GridInfo{T}, gridbox::GridBox{T}, chebcoefs::NTuple{3, ChebCoef{T}}) where{T}

    @assert length(qs) == length(poses)
    grid_revise_image!(gridbox)

    for i in 1:length(qs)
        interpolate_single_direct!(qs[i], poses[i], gridinfo, gridbox, chebcoefs)
    end

    return nothing
end