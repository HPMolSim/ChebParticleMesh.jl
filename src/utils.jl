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

