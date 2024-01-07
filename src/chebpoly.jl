"""
function ChebCoef(f::Function, h::T, w::Int, ν::Int) where{T<:Union{Float32, Float64}}

    This function creates a ChebCoef object for the Chebyshev approximation of function f.

    f: the function to be interpolated
    h: the grid spacing
    w: the number of grid points for the window function cutoff
    ν: the number of Chebyshev points
"""
function ChebCoef(f::Function, h::T, w::Int, ν::Int) where{T<:Union{Float32, Float64}}
    coef = pwcheb_coef(f, h, w, ν)
    return ChebCoef{T}(f, h, w, ν, coef)
end

function pwcheb_coef(f::Function, h::T, window_cut::Int, ν::Int) where{T<:Union{Float32, Float64}}
    window_width = T(window_cut * h)
    P = 2 * window_cut + 2
    C = zeros(T, P, ν)

    x = chebyshevpoints(T, ν, Val(1))

    Vdm = zeros(T, (ν, ν))
    for i in 1:ν
        for j in 1:ν
            Vdm[i, j] = x[i]^(j - 1)
        end
    end

    for i in 1:P
        if i == 1
            center = - window_width - h / 4
            scale = h / 4
        elseif i == P
            center = window_width + h / 4
            scale = h / 4
        else
            center = - window_width + (i - 1.5) * h
            scale = h / 2
        end

        b = f.(center .+ scale .* x)
        C[i, :] = Vdm \ b
    end

    return C
end

function pwcheb_eval!(x0::T, cheb_approx::Array{T, 1}, chebcoef::ChebCoef{T}) where{T<:Union{Float32, Float64}}

    h = chebcoef.h
    C = chebcoef.coef
    @assert (abs(x0) ≤ h / 2) | (abs(x0) ≈ h / 2)
    @assert size(cheb_approx, 1) == size(C, 1) - 1

    if x0 ≥ zero(T)
        cheb_approx[1] = horner((h / T(4) - x0) / (h / T(4)), @view C[1, :])
        for i in 2:size(C, 1) - 1
            cheb_approx[i] = horner((h / T(2) - x0) / (h / T(2)), @view C[i, :])
        end
    else
        cheb_approx[end] = horner(- (x0 + h / T(4)) / (h / T(4)), @view C[end, :])
        for i in 1:size(C, 1) - 2
            cheb_approx[i] = horner( - (x0 + h / T(2)) / (h / T(2)), @view C[i + 1, :])
        end
    end

    return cheb_approx
end

function pwcheb_eval(x0::T, chebcoef::ChebCoef{T}) where{T<:Union{Float32, Float64}}

    cheb_approx = zeros(T, size(chebcoef.coef, 1) - 1)
    cheb_approx = pwcheb_eval!(x0, cheb_approx, chebcoef)
    
    return cheb_approx
end

function f_eval!(x0::T, cheb_approx::Array{T, 1}, chebcoef::ChebCoef{T}) where{T<:Union{Float32, Float64}}

    h = chebcoef.h
    C = chebcoef.coef
    w = chebcoef.w
    @assert (abs(x0) ≤ h / 2) | (abs(x0) ≈ h / 2)
    @assert size(cheb_approx, 1) == size(C, 1) - 1

    for i in - w : w
        x = - x0 + i * h
        cheb_approx[i + w + 1] = chebcoef.f(x)
    end

    return cheb_approx
end

function f_eval(x0::T, chebcoef::ChebCoef{T}) where{T<:Union{Float32, Float64}}

    cheb_approx = zeros(T, 2 * chebcoef.w + 1)
    cheb_approx = f_eval!(x0, cheb_approx, chebcoef)
    
    return cheb_approx
end