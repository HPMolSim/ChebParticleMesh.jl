using ChebParticleMesh
using Test

@testset "ChebParticleMesh.jl" begin
    include("utils.jl")
    include("chebpoly.jl")
    include("index.jl")
    include("interpolate.jl")
    include("gather.jl")
end
