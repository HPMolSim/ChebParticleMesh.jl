@testset "pos/id transform" begin
    N_real = (20, 30, 40)
    w = (12, 23, 3)
    extra_pad = (24, 15, 6)
    L = (100.0, 100.0, 100.0)

    for periodicity_x in [true, false]
        for periodicity_y in [true, false]
            for periodicity_z in [true, false]
                periodicity = (periodicity_x, periodicity_y, periodicity_z)
                @testset "peridoicity is $(periodicity)" begin
                    gridinfo = GridInfo(N_real, w, periodicity, extra_pad, L)
                    for _ = 1:10
                        pos = Point(L[1] * rand(), L[2] * rand(), L[3] * rand())

                        id_image = image_grid_id(pos, gridinfo)
                        pos_trans_image = image_grid_pos(id_image, gridinfo)
                        for i in 1:3
                            @test abs(pos[i] - pos_trans_image[i]) ≤ gridinfo.h[i] / 2
                        end
                        
                        id_pad = pad_grid_id(pos, gridinfo)
                        pos_trans_pad = pad_grid_pos(id_pad, gridinfo)
                        for i in 1:3
                            @test abs(pos[i] - pos_trans_pad[i]) ≤ gridinfo.h[i] / 2
                        end
                    end
                end
            end
        end
    end
end

@testset "grid transform" begin
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