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

    for periodicity_x in [true, false], periodicity_y in [true, false], periodicity_z in [true, false]
        periodicity = (periodicity_x, periodicity_y, periodicity_z)
        @testset "peridoicity is $(periodicity)" begin
            gridinfo = GridInfo(N_real, w, periodicity, extra_pad, L)
            gb = GridBox(gridinfo)

            gb.image_grid .+= rand(gridinfo.N_image...)

            for i in 1:gridinfo.N_pad[1], j in 1:gridinfo.N_pad[2], k in 1:gridinfo.N_pad[3]
                id_image = id_pad2image(PadIndex((i, j, k)), gridinfo)
                for id_image_i in id_image
                    @test gb.pad_grid[i, j, k] == gb.image_grid[id_image_i.id...]
                end
            end

            for i in 1:gridinfo.N_image[1], j in 1:gridinfo.N_image[2], k in 1:gridinfo.N_image[3]
                id_pad = id_image2pad(ImageIndex((i, j, k)), gridinfo)
                @test gb.image_grid[i, j, k] == gb.pad_grid[id_pad.id...]
            end
        end
    end
end