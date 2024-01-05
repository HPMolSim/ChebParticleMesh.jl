@testset "interpolation by chebpoly" begin
    f = x -> exp(-x^2 * 10.0)
    h = 0.1
    w = 10
    ν = 10
    c = ChebCoef(f, h, w, ν)
    for x0 in [- h / 2 : 0.01 : h / 2...]
        p = [i * h for i in -w:w]
        @test pwcheb_eval(x0, c) ≈ f.(p .- x0)
    end
end