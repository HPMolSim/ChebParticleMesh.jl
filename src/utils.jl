"""
    function horner(x, p)
    
    Evaluate a polynomial using Horner's method.
"""
@inline function horner(x::T, p::AbstractArray{T, 1}) where{T}
    ex = p[end]
    for i = length(p)-1:-1:1
        ex = p[i] + x * ex
    end
    return ex
end

"""
function fft!(gridbox::GridBox{N, T}) where{N, T}

    calculate the inplace fft on the padded grid
"""
function FFTW.fft!(gridbox::GridBox{N, T}) where{N, T}

    fft!(gridbox.pad_grid)

    return gridbox
end

"""
function ifft!(gridbox::GridBox{N, T}) where{N, T}

    calculate the inplace ifft on the padded grid
"""
function FFTW.ifft!(gridbox::GridBox{N, T}) where{N, T}

    ifft!(gridbox.pad_grid)

    return gridbox
end

function funcpack(f::Function, args::Vector)
    g = x -> f(x, args...)
    return g
end

"""
    function Wkb(x::T, width::T, β::T) where{T<:Real}

    WKB kernel function, used as window function for the Chebyshev interpolation
"""
function Wkb(x::T, width::T, β::T) where{T<:Real}
    return T(besseli(0, β * sqrt(abs(one(T) - (x / width)^2))) / besseli(0, β) * (abs(x) <= width))
end

"""
    function FWkb(k::T, width::T, β::T) where{T<:Real}

    Fourier transform of Wkb kernel function
"""
function FWkb(k::T, width::T, β::T) where{T}
    return T(2 * width * sinh(sqrt(β^2 - (k * width)^2)) / (besseli(0, β) * sqrt(β^2 - (k * width)^2)))
end