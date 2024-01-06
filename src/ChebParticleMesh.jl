module ChebParticleMesh

using LinearAlgebra, ExTinyMD, FastTransforms, SpecialFunctions, LoopVectorization

export horner
export ChebCoef, pwcheb_eval, pwcheb_eval!
export GridInfo, PadIndex, ImageIndex, id_image2pad, id_pad2image, grid_revise_image!, grid_revise_pad!, grid_image2pad!, grid_pad2image!

include("types.jl")
include("utils.jl")

include("chebpoly.jl")
include("grid.jl")
include("index.jl")

include("interpolation.jl")
include("scaling.jl")
include("gather.jl")

end
