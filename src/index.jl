function nearest_grid_id(x::NTuple{N, T}, h::NTuple{N, T}) where{N, T}

    id = Int.(round.((x .- h./2) ./ h)) .+ 1

    return id
end

function pad_grid_id(x::NTuple{N, T}, gridinfo::GridInfo{N, T}) where{N, T}
    
    id = nearest_grid_id(x, gridinfo.h)

    return PadIndex(id .+ gridinfo.pad)
end

function pad_grid_pos(id::PadIndex{N}, gridinfo::GridInfo{N, T}) where{N, T}

    pos = (id.id .- gridinfo.pad .- T(0.5)) .* gridinfo.h

    return pos
end

function image_grid_id(x::NTuple{N, T}, gridinfo::GridInfo{N, T}) where{N, T}
    
    id = nearest_grid_id(x, gridinfo.h)

    return ImageIndex(id .+ gridinfo.image)
end

function image_grid_pos(id::ImageIndex{N}, gridinfo::GridInfo{N, T}) where{N, T}

    pos = (id.id .- gridinfo.image .- T(0.5)) .* gridinfo.h

    return pos
end

function id_image2pad_single(id_image_i::Int, N_i::T, image_i::Int, pad_i::Int, periodicity::Bool) where{T}

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

function id_image2pad(image::ImageIndex{N}, gridinfo::GridInfo{N, T}) where{N, T}

    pad_id = tuple([id_image2pad_single(image.id[i], gridinfo.N_real[i], gridinfo.image[i], gridinfo.pad[i], gridinfo.periodicity[i]) for i in 1:N]...)

    # px = id_image2pad_single(ix, gridinfo.N_real[1], gridinfo.image[1], gridinfo.pad[1], gridinfo.periodicity[1])
    # py = id_image2pad_single(iy, gridinfo.N_real[2], gridinfo.image[2], gridinfo.pad[2], gridinfo.periodicity[2])
    # pz = id_image2pad_single(iz, gridinfo.N_real[3], gridinfo.image[3], gridinfo.pad[3], gridinfo.periodicity[3])

    return PadIndex(pad_id)
end

function id_pad2image_single(id_pad_i::Int, N_i::T, image_i::Int, pad_i::Int, periodicity::Bool) where{T}

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

function id_pad2image(pad::PadIndex{N}, gridinfo::GridInfo{N, T}) where{N, T}

    image_id = [id_pad2image_single(pad.id[i], gridinfo.N_real[i], gridinfo.image[i], gridinfo.pad[i], gridinfo.periodicity[i]) for i in 1:N]

    all_image_id = Vector{ImageIndex{N}}()
    
    for i in Iterators.product(image_id...)
        push!(all_image_id, ImageIndex(i))
    end

    return all_image_id
end