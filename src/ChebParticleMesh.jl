module ChebParticleMesh

using LinearAlgebra, FastTransforms, SpecialFunctions, LoopVectorization, FFTW

export horner
export ChebCoef, pwcheb_eval, pwcheb_eval!, f_eval, f_eval!
export GridInfo, GridBox, PadIndex, ImageIndex
export id_image2pad, id_pad2image, grid_revise_pad!
export nearest_grid_id, pad_grid_id, pad_grid_pos, image_grid_id, image_grid_pos
export interpolate_single!, interpolate!, interpolate_single_direct!, interpolate_direct!

include("types.jl")
include("utils.jl")

include("chebpoly.jl")
include("grid.jl")
include("index.jl")

include("interpolate.jl")
include("scaling.jl")
include("gather.jl")

end
