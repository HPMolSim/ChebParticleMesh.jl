@testset "gather 1d" begin

    q = 0.8

    for pos in [tuple(0.01), tuple(0.321), tuple(0.5)]
        for periodicity in [true, false]
            @testset "pos is $(pos), peridoicity is $(periodicity)" begin
                gi = GridInfo(tuple(100), tuple(10), tuple(periodicity), tuple(10), tuple(10.0))
                f_window = [x -> Wkb(x, (gi.w[i] + 0.5) * gi.h[i], 5.0 * gi.w[i]) for i in 1:1]
                cheb_coefs = tuple([ChebCoef(f_window[i], gi.h[i], gi.w[i], 10) for i in 1:1]...)

                gb = GridBox(gi)
                gb.pad_grid .+= rand(eltype(gb.pad_grid), size(gb.pad_grid))

                energy_cheb = gather_single(q, pos, gi, gb, cheb_coefs)
                energy_direct = gather_single_direct(q, pos, gi, gb, cheb_coefs)
                energy_naive = gather_single_naive(q, pos, gi, gb, cheb_coefs)

                @test energy_cheb ≈ energy_naive ≈ energy_direct
            end
        end
    end
end

@testset "gather 2d" begin

    q = 0.8

    for pos in [tuple(0.01, 0.01), tuple(0.321, 5.232), tuple(0.5, 2.0)]
        for periodicity in [(true, true), (true, false), (false, true), (false, false)]
            @testset "pos is $(pos), peridoicity is $(periodicity)" begin
                gi = GridInfo(tuple(100, 110), tuple(10, 8), periodicity, (10, 10), tuple(9.0, 11.0))
                f_window = [x -> Wkb(x, (gi.w[i] + 0.5) * gi.h[i], 5.0 * gi.w[i]) for i in 1:2]
                cheb_coefs = tuple([ChebCoef(f_window[i], gi.h[i], gi.w[i], 10) for i in 1:2]...)

                gb = GridBox(gi)
                gb.pad_grid .+= rand(eltype(gb.pad_grid), size(gb.pad_grid))

                energy_cheb = gather_single(q, pos, gi, gb, cheb_coefs)
                energy_direct = gather_single_direct(q, pos, gi, gb, cheb_coefs)
                energy_naive = gather_single_naive(q, pos, gi, gb, cheb_coefs)

                @test energy_cheb ≈ energy_naive ≈ energy_direct
            end
        end
    end

end