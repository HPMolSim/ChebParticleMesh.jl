@testset "grid transform: image to pad" begin
    N_real = (5, 6, 7)
    w = (1, 3, 4)
    extra_pad = (4, 5, 6)
    L = (100.0, 100.0, 100.0)

    for periodicity_x in [true, false]
        for periodicity_y in [true, false]
            for periodicity_z in [true, false]
                periodicity = (periodicity_x, periodicity_y, periodicity_z)
                @testset "peridoicity is $(periodicity)" begin
                    gridinfo = GridInfo(N_real, w, periodicity, extra_pad, L)
                    for transinfo_i in gridinfo.trans_info
                        ix, iy, iz = transinfo_i[1]
                        px, py, pz = transinfo_i[2]
                        sx, sy, sz = transinfo_i[3]
                        for i in 1:sx
                            for j in 1:sy
                                for k in 1:sz
                                    @test id_image2pad(ImageIndex(ix + i, iy + j, iz + k), gridinfo) == PadIndex(px + i, py + j, pz + k)
                                    @test ImageIndex(ix + i, iy + j, iz + k) ∈ id_pad2image(PadIndex(px + i, py + j, pz + k), gridinfo)
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end