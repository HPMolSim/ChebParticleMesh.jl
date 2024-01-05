struct ChebCoef{T}
    h::T
    w::Int
    ν::Int
    coef::Array{T, 2}
end

struct GridInfo{T}
    N_real::NTuple{3, Int}
    w::NTuple{3, Int}
    periodicity::NTuple{3, Bool}
    image::NTuple{3, Int}
    pad::NTuple{3, Int}
    L::NTuple{3, T}
    h::NTuple{3, T}

    N_image::NTuple{3, Int}
    N_pad::NTuple{3, Int}
    
    trans_info::Vector{Tuple{NTuple{3, Int}, NTuple{3, Int}, NTuple{3, Int}}}
end

mutable struct GridBox{T}
    pad_grid::Array{T, 3}
    image_grid::Array{T, 3}
end

abstract type AbstractIndex end

struct PadIndex <: AbstractIndex
    idx::Int
    idy::Int
    idz::Int
end

struct ImageIndex <: AbstractIndex
    idx::Int
    idy::Int
    idz::Int
end