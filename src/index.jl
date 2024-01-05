"""

"""
function nearest_grid_id(x::Point{3, T}, h::NTuple{3, T}) where{T}

    idx = Int(round((x[1] - h[1] / 2) / h[1])) + 1
    idy = Int(round((x[2] - h[2] / 2) / h[2])) + 1
    idz = Int(round((x[3] - h[3] / 2) / h[3])) + 1

    return idx, idy, idz
end

function pad_grid_id(x::Point{3, T}, gridinfo::GridInfo{T}) where{T}
    
    idx, idy, idz = nearest_grid_id(x, gridinfo.h)

    idx = idx + gridinfo.pad[1]
    idy = idy + gridinfo.pad[2]
    idz = idz + gridinfo.pad[3]

    return PadIndex(idx, idy, idz)
end

function image_grid_id(x::Point{3, T}, gridinfo::GridInfo{T}) where{T}
    
    idx, idy, idz = nearest_grid_id(x, gridinfo.h)

    idx = idx + gridinfo.w[1] + 1
    idy = idy + gridinfo.w[2] + 1
    idz = idz + gridinfo.w[3] + 1

    return ImageIndex(idx, idy, idz)
end

function image2pad_single(id_image_i::Int, N_i::T, image_i::Int, pad_i::Int, periodicity::Bool) where{T}

    id_pad_i = 0
    if id_image_i ≤ image_i
        if periodicity == true
            id_pad_i = id_image_i - image_i + N_i
        else
            id_pad_i = id_image_i + pad_i - image_i
        end
    elseif id_image_i ≥ N_i + image_i + 1
        if periodicity == true
            id_pad_i = id_image_i - image_i - N_i
        else
            id_pad_i = id_image_i + pad_i - image_i
        end
    else
        id_pad_i = id_image_i + pad_i - image_i
    end

    return id_pad_i
end

function image2pad(image::ImageIndex, gridinfo::GridInfo{T}) where{T}

    ix, iy, iz = image.idx, image.idy, image.idz

    px = image2pad_single(ix, gridinfo.N_real[1], gridinfo.image[1], gridinfo.pad[1], gridinfo.periodicity[1])
    py = image2pad_single(iy, gridinfo.N_real[2], gridinfo.image[2], gridinfo.pad[2], gridinfo.periodicity[2])
    pz = image2pad_single(iz, gridinfo.N_real[3], gridinfo.image[3], gridinfo.pad[3], gridinfo.periodicity[3])

    return PadIndex(px, py, pz)
end

function pad2image_single(id_pad_i::Int, N_i::T, image_i::Int, pad_i::Int, periodicity::Bool) where{T}

    if periodicity == false
        return [id_pad_i - pad_i + image_i]
    else
        if (id_pad_i ≤ image_i) & (id_pad_i ≥ N_i - image_i + 1)
            return [id_pad_i + image_i, id_pad_i + N_i + image_i, id_pad_i - N_i + image_i]
        elseif id_pad_i ≤ image_i
            return [id_pad_i + image_i, id_pad_i + N_i + image_i]
        elseif id_pad_i ≥ N_i - image_i + 1
            return [id_pad_i + image_i, id_pad_i - N_i + image_i]
        else
            return [id_pad_i + image_i]
        end
    end
end

function pad2image(pad::PadIndex, gridinfo::GridInfo{T}) where{T}

    px, py, pz = pad.idx, pad.idy, pad.idz

    ixs = pad2image_single(px, gridinfo.N_real[1], gridinfo.image[1], gridinfo.pad[1], gridinfo.periodicity[1])
    iys = pad2image_single(py, gridinfo.N_real[2], gridinfo.image[2], gridinfo.pad[2], gridinfo.periodicity[2])
    izs = pad2image_single(pz, gridinfo.N_real[3], gridinfo.image[3], gridinfo.pad[3], gridinfo.periodicity[3])

    return [ImageIndex(ix, iy, iz) for ix in ixs, iy in iys, iz in izs]
end