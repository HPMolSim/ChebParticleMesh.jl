@testset "pos/id transform 2D" begin
    N_real = (20, 30)
    w = (12, 23)
    extra_pad = (24, 15)
    L = (100.0, 100.0)

    for periodicity_x in [true, false]
        for periodicity_y in [true, false]
            periodicity = (periodicity_x, periodicity_y)
            @testset "peridoicity is $(periodicity)" begin
                gridinfo = GridInfo(N_real, w, periodicity, extra_pad, L)
                for _ = 1:10
                    pos = (L[1] * rand(), L[2] * rand())

                    id_image = image_grid_id(pos, gridinfo)
                    pos_trans_image = image_grid_pos(id_image, gridinfo)
                    for i in 1:2
                        @test abs(pos[i] - pos_trans_image[i]) ≤ gridinfo.h[i] / 2
                    end
                    
                    id_pad = pad_grid_id(pos, gridinfo)
                    pos_trans_pad = pad_grid_pos(id_pad, gridinfo)
                    for i in 1:2
                        @test abs(pos[i] - pos_trans_pad[i]) ≤ gridinfo.h[i] / 2
                    end
                end
            end
        end
    end
end

@testset "pos/id transform 3D" begin
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
                        pos = (L[1] * rand(), L[2] * rand(), L[3] * rand())

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
                        image_region, pad_region = transinfo_i
                        for (ix, px) in zip(image_region[1], pad_region[1])
                            for (iy, py) in zip(image_region[2], pad_region[2])
                                for (iz, pz) in zip(image_region[3], pad_region[3])
                                    @test id_image2pad(ImageIndex((ix, iy, iz)), gridinfo) == PadIndex((px, py, pz))
                                    @test ImageIndex((ix, iy, iz)) ∈ id_pad2image(PadIndex((px, py, pz)), gridinfo)
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end