module ChebParticleMesh

using LinearAlgebra, ExTinyMD, FastTransforms, SpecialFunctions

export horner
export ChebCoef, pwcheb_eval, pwcheb_eval!
export GridInfo, PadIndex, ImageIndex, image2pad, pad2image

include("types.jl")
include("utils.jl")

include("chebpoly.jl")
include("grid.jl")
include("index.jl")

end
