function ScalingFactor(f::Function, gridinfo::GridInfo{N, T}) where{N, T}
    factors = zeros(Complex{T}, gridinfo.N_pad...)

    for i in Iterators.product([1:gridinfo.N_pad[d] for d in 1:N]...)
        k = [gridinfo.k[d][i[d]] for d in 1:N]
        factors[i...] = Complex{T}(f(k...))
    end

    return ScalingFactor{N, T}(f, factors)
end

function scale!(gridbox::GridBox{N, T}, scalingfactor::ScalingFactor{N, T}) where{N, T}

    gridbox.pad_grid .*= scalingfactor.factors

    return gridbox
end