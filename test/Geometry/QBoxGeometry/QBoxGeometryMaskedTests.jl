module QBoxGeometryMaskedTests
using Mantis
using Test

# 1D base geometry, 2 elements
cg = Geometry.CartesianGeometry((LinRange(0, 1, 3),)) 
# evaluation mask with 2 base elements
E = Geometry.trivial_evaluation_mask(1, 2)
# Masked geometry
mg = Geometry.MaskedGeometry(cg, E)
qbox_size1d = (2,)
num_sub1d = (3,)
qbg_mp1d = Geometry.QBoxGeometry_refine(mg, qbox_size1d, num_sub1d)
#Plot.plot(mg; vtk_filename="MaskGeo 1D")

# # 2D base geometry
# cg2 = Geometry.CartesianGeometry((
#     LinRange(0, 1, 3),   # 2 elementen in x
#     LinRange(0, 1, 4),   # 3 elementen in y
# ))
# # evaluation mask with 2 base elements
# E2 = Geometry.trivial_evaluation_mask(2, 2)
# # Masked geometry
# mg2 = Geometry.MaskedGeometry(cg2, E2)
# qbox_size2d=(2,1)
# num_sub2d=(2,2)
# qbg_mp2d = Geometry.QBoxGeometry_from_existing(mg2, qbox_size2d, num_sub2d)
#Plot.plot(mg2; vtk_filename="MaskGeo 2D")

@testset "QBoxGeometry Masked geom 1d" begin
    qid, lvl, pid = Geometry.get_qbox_id_hier(1, qbg_mp1d)
    @test lvl == 1
    @test qid == 1
    @test pid == 1
    
    qid, lvl, pid = Geometry.get_qbox_id_hier(3, qbg_mp1d)
    @test lvl == 1
    @test qid == 2
    @test pid == 1

    @test Geometry.get_qbox_element_ids(1, 1, qbg_mp1d, 2) == [3,4]

    Geometry.refine_qbox!(qbg_mp1d, 1, 1, 1)

    children_lvl2 = Geometry.get_child_qbox_ids(1, 1, qbg_mp1d, 2)
    @test children_lvl2 == [4,5,6]

    # hier_id = 1 is now the first active element on level 1
    qid, lvl, pid = Geometry.get_qbox_id_hier(1, qbg_mp1d)
    @test lvl == 1
    @test qid == 2
    @test pid == 1

    qid, lvl, pid = Geometry.get_qbox_id_hier(6, qbg_mp1d)
    @test lvl == 2
    @test qid == 2
    @test pid == 1


    @test Geometry.get_qbox_element_ids(2, 1, qbg_mp1d, 3) == [5,6]
end

end