struct ChebCoef{T}
    f::Function
    h::T
    w::Int
    ν::Int
    coef::Array{T, 2}
end

struct GridInfo{N, T}
    N_real::NTuple{N, Int}
    w::NTuple{N, Int}
    periodicity::NTuple{N, Bool}
    image::NTuple{N, Int}
    pad::NTuple{N, Int}
    L::NTuple{N, T}
    h::NTuple{N, T}

    N_image::NTuple{N, Int}
    N_pad::NTuple{N, Int}
    index_list::Vector{Vector{Int}}
    iter_list::Iterators.ProductIterator{NTuple{N, UnitRange{Int64}}}

    k::Vector{Array{T, 1}}
end

mutable struct GridBox{N, T}
    pad_grid::Array{Complex{T}, N}
    image_grid::SubArray{Complex{T}, N}

    cheb_value::Vector{Array{T, 1}}
end

abstract type AbstractIndex end

struct PadIndex{N} <: AbstractIndex
    id::NTuple{N, Int}
end

struct ImageIndex{N} <: AbstractIndex
    id::NTuple{N, Int}
end

struct ScalingFactor{N, T}
    f::Function
    factors::Array{Complex{T}, N}
end