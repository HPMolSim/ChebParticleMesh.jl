using ChebParticleMesh
using Test
using ExTinyMD

@testset "ChebParticleMesh.jl" begin
    include("utils.jl")
    include("chebpoly.jl")
    include("index.jl")
    include("interpolate.jl")
end
