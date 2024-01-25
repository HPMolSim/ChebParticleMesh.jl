@testset "interpolate 1d" begin

    q = 0.8

    for pos in [tuple(0.01), tuple(0.321), tuple(0.5)]
        for periodicity in [true, false]
            @testset "pos is $(pos), peridoicity is $(periodicity)" begin
                gi = GridInfo(tuple(100), tuple(10), tuple(periodicity), tuple(10), tuple(10.0))
                f_window = [x -> Wkb(x, (gi.w[i] + 0.5) * gi.h[i], 5.0 * gi.w[i]) for i in 1:1]
                cheb_coefs = tuple([ChebCoef(f_window[i], gi.h[i], gi.w[i], 10) for i in 1:1]...)

                gb_cheb = GridBox(gi)
                gb_direct = GridBox(gi)
                gb_naive = GridBox(gi)

                interpolate_single!(q, pos, gi, gb_cheb, cheb_coefs)
                interpolate_single_direct!(q, pos, gi, gb_direct, cheb_coefs)
                interpolate_single_naive!(q, pos, gi, gb_naive, cheb_coefs)
            
                @test gb_cheb.pad_grid ≈ gb_naive.pad_grid
                @test gb_cheb.pad_grid ≈ gb_direct.pad_grid
                @test gb_cheb.image_grid ≈ gb_naive.image_grid
                @test gb_cheb.image_grid ≈ gb_direct.image_grid
            end
        end
    end

end

@testset "interpolate 2d" begin

    q = 0.8

    for pos in [tuple(0.01, 0.01), tuple(0.321, 5.232), tuple(0.5, 2.0)]
        for periodicity in [(true, true), (true, false), (false, true), (false, false)]
            @testset "pos is $(pos), peridoicity is $(periodicity)" begin
                gi = GridInfo(tuple(100, 110), tuple(10, 8), periodicity, (10, 10), tuple(9.0, 11.0))
                f_window = [x -> Wkb(x, (gi.w[i] + 0.5) * gi.h[i], 5.0 * gi.w[i]) for i in 1:2]
                cheb_coefs = tuple([ChebCoef(f_window[i], gi.h[i], gi.w[i], 10) for i in 1:2]...)

                gb_cheb = GridBox(gi)
                gb_direct = GridBox(gi)
                gb_naive = GridBox(gi)

                interpolate_single!(q, pos, gi, gb_cheb, cheb_coefs)
                interpolate_single_direct!(q, pos, gi, gb_direct, cheb_coefs)
                interpolate_single_naive!(q, pos, gi, gb_naive, cheb_coefs)
            
                @test gb_cheb.pad_grid ≈ gb_naive.pad_grid
                @test gb_cheb.pad_grid ≈ gb_direct.pad_grid
                @test gb_cheb.image_grid ≈ gb_naive.image_grid
                @test gb_cheb.image_grid ≈ gb_direct.image_grid
            end
        end
    end

end