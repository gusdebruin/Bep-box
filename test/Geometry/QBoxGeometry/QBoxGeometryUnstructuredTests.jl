module QBoxGeometryUnstructuredTests
using Mantis
using Test

############################################################################################
#                                          Setup                                           #
############################################################################################

# Test UnstructuredGeometry ------------------------------------------------

# Reduction test, single-patch, single element, 1D.
geometry_1 = Geometry.UnstructuredGeometry((Geometry.CartesianGeometry(([-1, 1],)),))
qbox_size_1 = (4,)
num_sub_1 = (2,)
qbp_1 = Geometry.QBoxGeometry_refine(geometry_1, qbox_size_1, num_sub_1)


# LinRange input. Single-patch, 2D.
cg1 = Geometry.CartesianGeometry((LinRange(0.5, 2.5, 5), LinRange(-0.75, 0.75, 3)))
cg2 = Geometry.CartesianGeometry((LinRange(2.5, 5.0, 6), LinRange(-0.75, 0.75, 5)))
geometry_2 = Geometry.UnstructuredGeometry((cg1, cg2))
qbox_size_2 = (1,2)
num_sub_2 = (2,1)
qbp_2 = Geometry.QBoxGeometry_from_existing(geometry_2, qbox_size_2, num_sub_2)



############################################################################################
#                                             Tests                                        #
############################################################################################

@testset "QBoxGeometry Unstructured geom, two patches" begin
     # element 1 → qbox 1
    qid, lvl, pid = Geometry.get_qbox_id_hier(1, qbp_2)
    @test lvl == 1
    @test qid == 1
    @test pid == 1
    
    qid, lvl, pid = Geometry.get_qbox_id_hier(6, qbp_2)
    @test lvl == 1
    @test qid == 2
    @test pid == 1

    # element 28  → qbox 4
    qid, lvl, pid = Geometry.get_qbox_id_hier(22, qbp_2)
    @test lvl == 1
    @test qid == 9
    @test pid == 2

    # qbox 4 of patch 2 contains 4 elements on level 1 (not that this returns level ids)
    @test Geometry.get_qbox_element_ids(1, 2, qbp_2, 3) == [11,16]

    # level 2:
    Geometry.refine_qbox!(qbp_2, 1, 1, 1)
    Geometry.refine_qbox!(qbp_2, 1, 1, 2)
    Geometry.refine_qbox!(qbp_2, 1, 2, 1)
    Geometry.refine_qbox!(qbp_2, 1, 2, 6)

    # level 3: 
    Geometry.refine_qbox!(qbp_2, 2, 1, 2)
    Geometry.refine_qbox!(qbp_2, 2, 2, 2)
    
    #Plot.plot(qbp_2; vtk_filename="After QBOX refinement Unstructured, multiple patches")

    children_lvl2 = Geometry.get_child_qbox_ids(1, 1, qbp_2, 4)
    @test children_lvl2 == [7,8]

    children_lvl2 = Geometry.get_child_qbox_ids(2, 2, qbp_2, 11)
    @test children_lvl2 == [21,22]

    # hier_id = 1 is now the first active element on level 1
    qid, lvl, pid = Geometry.get_qbox_id_hier(1, qbp_2)
    @test lvl == 1
    @test qid == 3
    @test pid == 1

    qid, lvl, pid = Geometry.get_qbox_id_hier(31, qbp_2)
    @test lvl == 2
    @test qid == 12
    @test pid == 2

    qid, lvl, pid = Geometry.get_qbox_id_hier(24, qbp_2)
    @test lvl == 2
    @test qid == 3
    @test pid == 1

    qid, lvl, pid = Geometry.get_qbox_id_hier(22, qbp_2)
    @test lvl == 2
    @test qid == 1
    @test pid == 1

    qid, lvl, pid = Geometry.get_qbox_id_hier(8, qbp_2)
    @test lvl == 1
    @test qid == 5
    @test pid == 2

    qid, lvl, pid = Geometry.get_qbox_id_hier(40, qbp_2)
    @test lvl == 3
    @test qid == 4
    @test pid == 2

    @test Geometry.get_qbox_element_ids(3, 2, qbp_2, 3) == [35,55]

end

@testset "QBoxGeometry Unstructured geom, one element" begin
     # element 1 → qbox 1
    qid, lvl, pid = Geometry.get_qbox_id_hier(1, qbp_1)
    @test lvl == 1
    @test qid == 1
    @test pid == 1
    
    qid, lvl, pid = Geometry.get_qbox_id_hier(4, qbp_1)
    @test lvl == 1
    @test qid == 1
    @test pid == 1


    # qbox 4 of patch 2 contains 4 elements on level 1 (not that this returns level ids)
    @test Geometry.get_qbox_element_ids(1, 1, qbp_1, 1) == [1,2,3,4]

    # level 2:
    Geometry.refine_qbox!(qbp_1, 1, 1, 1)

    # level 3: 
    Geometry.refine_qbox!(qbp_1, 2, 1, 1)
    
    children_lvl2 = Geometry.get_child_qbox_ids(1, 1, qbp_1, 1)
    @test children_lvl2 == [1,2]

    children_lvl2 = Geometry.get_child_qbox_ids(2, 1, qbp_1, 1)
    @test children_lvl2 == [1,2]

    # hier_id = 1 is now the first active element on level 1
    qid, lvl, pid = Geometry.get_qbox_id_hier(1, qbp_1)
    @test lvl == 2
    @test qid == 2
    @test pid == 1

    qid, lvl, pid = Geometry.get_qbox_id_hier(9, qbp_1)
    @test lvl == 3
    @test qid == 2
    @test pid == 1

    @test Geometry.get_qbox_element_ids(3, 1, qbp_1, 1) == [1,2,3,4]

end

end
