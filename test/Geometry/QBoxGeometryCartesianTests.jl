module QBoxGeometryTests

using Mantis
using Test

############################################################################################
#                                          Setup                                           #
############################################################################################

# 1D QBox two-level hierarchy
# create a 1D geometry with 4 elements : [[0,1],[1,2],[2,3],[3,4]] for level 1
geom_lvl1_1D = Geometry.CartesianGeometry((0.0:2.0:4.0,))
n_elements_lvl1_1D = Geometry.get_num_elements(geom_lvl1_1D)
qbox_size_1D = (2,)
num_subdivisions_1D = (2,)
qbg_1D = Geometry.QBoxGeometry_refine(geom_lvl1_1D, qbox_size_1D, num_subdivisions_1D)

#2D QBox three-level hierarchy
geom_lvl1_2D = Geometry.CartesianGeometry((
    0.0:2.0:4.0,
    0.0:2.0:4.0
))
n1_2D = Geometry.get_num_elements(geom_lvl1_2D)   # = 16
qbox_size_2D = (2,2)
num_subdivisions_2D = (2,2)
qbg_2D = Geometry.QBoxGeometry_refine(geom_lvl1_2D, qbox_size_2D, num_subdivisions_2D)

# 2D multi patch geometry
geom_lvl1_2D_mp = Geometry.CartesianGeometry(((LinRange(0.5, 2.5, 3), LinRange(-0.75, 0.75, 2)),(LinRange(2.5, 5.5, 4), LinRange(-0.75, 0.75, 3))))
qbox_size_2D_mp = (2,2)
num_subdivisions_2D_mp = (2,2)
qbg_2D_mp = Geometry.QBoxGeometry_refine(geom_lvl1_2D_mp, qbox_size_2D_mp, num_subdivisions_2D_mp)

# 2D multi patch geometry from exisiting
geom_lvl1_2D_mp_e =Geometry.CartesianGeometry(((LinRange(0, 6, 7), LinRange(0, 4, 5)),(LinRange(6, 15, 10), LinRange(0, 4, 3))))
qbox_size_2D_mp_e = (3,2)
num_subdivisions_2D_mp_e = (2,3)
qbg_2D_mp_e = Geometry.QBoxGeometry_from_existing(geom_lvl1_2D_mp_e, qbox_size_2D_mp_e, num_subdivisions_2D_mp_e)
#Plot.plot(qbg_2D_mp_e; vtk_filename="test 2D multi-patch_e before QBox")

# 3D geometry
geom_3D_mp = Geometry.CartesianGeometry((
    (
        LinRange(0, 2, 3),   # x: 2 elements
        LinRange(0, 2, 3),   # y: 2 elements
        LinRange(0, 2, 3)    # z: 2 elements
    ),
    (
        LinRange(2, 5, 4),   # x: 3 elements
        LinRange(0, 2, 2),   # y: 1 element
        LinRange(0, 2, 3)    # z: 2 elements
    )
))
qbox_size_3D = (1,1,1)
num_subdivisions_3D = (2,2,2)
qbg_3D_mp = Geometry.QBoxGeometry_from_existing(geom_3D_mp, qbox_size_3D, num_subdivisions_3D)
#Plot.plot(geom_3D_mp; vtk_filename="Basic 3D multi-patch_e before QBox")
# testing of beginning with a QBoxGeometry or HierarchicalGeometry gives an error (they all gave the correct error)
# geom_lvl1_hier = Geometry.CartesianGeometry((
#     0.0:2.0:4.0,
#     0.0:2.0:4.0
# ))
# n1_2D_hier = Geometry.get_num_elements(geom_lvl1_hier)   # = 16
# qbox_size_2D_hier = (2,2)
# num_subdivisions_2D_hier = (2,2)

# geom_lvl2_hier = Geometry.CartesianGeometry((
#     0.0:1.0:4.0,
#     0.0:1.0:4.0
# ))

# # ActiveInfo: all elements of level 1 are active
# active_2D_hier = Hierarchy.ActiveInfo([collect(1:n1_2D_hier), Int[]])

# hier_g_test = Geometry.HierarchicalGeometry(
#     (geom_lvl1_hier, geom_lvl2_hier),
#     active_2D_hier
# )
# qbg_2D_test = Geometry.QBoxGeometry(hier_g_test, qbox_size_2D_hier, num_subdivisions_2D_hier)

# #test_hier_r = Geometry.QBoxGeometry_refine(hier_g_test, qbox_size_2D_hier, num_subdivisions_2D_hier)
# #test_qbg_r = Geometry.QBoxGeometry_refine(qbg_2D_test, qbox_size_2D_hier, num_subdivisions_2D_hier)
# #test_hier_e = Geometry.QBoxGeometry_from_existing(hier_g_test, qbox_size_2D_hier, num_subdivisions_2D_hier)
# #test_qbg_e = Geometry.QBoxGeometry_from_existing(qbg_2D_test, qbox_size_2D_hier, num_subdivisions_2D_hier)

############################################################################################
#                                       Basic Tests                                        #
############################################################################################

@testset "QBoxGeometry basic tests 1D refine" begin
    # element 1 → qbox 1
    qid, lvl, pid = Geometry.get_qbox_id_hier(1, qbg_1D)
    @test qid == 1
    @test lvl == 1
    @test pid == 1

    qid, lvl, pid = Geometry.get_qbox_id_hier(3, qbg_1D)
    @test qid == 2
    @test lvl == 1
    @test pid == 1

    # qbox 1 contains elements 1 and 2
    @test Geometry.get_qbox_element_ids(1, 1, qbg_1D, 1) == [1,2]

    active = Geometry.get_active_elements(Geometry.get_hierarchical_geometry(qbg_1D))
    @show Hierarchy.get_level_ids(active)
    # refine qbox 1
    Geometry.refine_qbox!(qbg_1D, 1, 1, 1)
    active = Geometry.get_active_elements(Geometry.get_hierarchical_geometry(qbg_1D))
    @show Hierarchy.get_level_ids(active)

    # children of qbox 1 (level 1 → level 2)
    children = Geometry.get_child_qbox_ids(1, 1, qbg_1D, 1)
    @test length(children) == 2
    @test children == [1,2]

    qid, lvl, pid = Geometry.get_qbox_id_hier(1, qbg_1D)
    @test qid == 2
    @test lvl == 1
    @test pid == 1

    qid, lvl, pid = Geometry.get_qbox_id_hier(3, qbg_1D)
    @test qid == 1
    @test lvl == 2
    @test pid == 1
end

@testset "QBoxGeometry basic tests 2D refine" begin
    # element 1 → qbox 1
    qid, lvl, pid = Geometry.get_qbox_id_hier(1, qbg_2D)
    @test lvl == 1
    @test qid == 1

    # element 6 (row 2, col 2) → qbox 1
    qid, lvl, pid = Geometry.get_qbox_id_hier(6, qbg_2D)
    @test lvl == 1
    @test qid == 1

    # element 10 (row 3, col 2) → qbox 3
    qid, lvl, pid = Geometry.get_qbox_id_hier(10, qbg_2D)
    @test lvl == 1
    @test qid == 3

    # qbox 1 contains 4 elements on level 1
    @test Geometry.get_qbox_element_ids(1, 1, qbg_2D, 1) == [1,2,5,6]


    Geometry.refine_qbox!(qbg_2D, 1, 1, 1)
    # Now active elements:
    # level 1: elements not in qbox 1
    # level 2: children of qbox 1

    # children of qbox 1 on level 2
    children_lvl2 = Geometry.get_child_qbox_ids(1, 1, qbg_2D, 2)
    @test length(children_lvl2) == 4
    @test children_lvl2 == [3,4,7,8]


    # hier_id = 1 is now the first active element on level 1
    qid, lvl, pid = Geometry.get_qbox_id_hier(1, qbg_2D)
    @test lvl == 1
    @test qid == 2

    qid, lvl, pid = Geometry.get_qbox_id_hier(12, qbg_2D)
    @test lvl == 1
    # hier_id = 5 should now be on level 2
    qid, lvl, pid = Geometry.get_qbox_id_hier(13, qbg_2D)
    @test lvl == 2
end

@testset "QBoxGeometry basic tests 2D Multi-patch refine" begin
    # element 1 → qbox 1
    qid, lvl, pid = Geometry.get_qbox_id_hier(1, qbg_2D_mp)
    @test lvl == 1
    @test qid == 1
    @test pid == 1
    

    # element 7 (row 2, col 3) → qbox 2
    qid, lvl, pid = Geometry.get_qbox_id_hier(7, qbg_2D_mp)
    @test lvl == 1
    @test qid == 2
    @test pid == 1

    # element 28  → qbox 4
    qid, lvl, pid = Geometry.get_qbox_id_hier(28, qbg_2D_mp)
    @test lvl == 1
    @test qid == 4
    @test pid == 2

    # qbox 4 of patch 2 contains 4 elements on level 1 (not that this returns level ids)
    @test Geometry.get_qbox_element_ids(1, 2, qbg_2D_mp, 4) == [21,22,27,28]

    Geometry.refine_qbox!(qbg_2D_mp, 1, 1, 1)
    Geometry.refine_qbox!(qbg_2D_mp, 1, 2, 2)
    # Now active elements:
    # level 1: elements not in qbox 1 of patch 1 and qbox 2 of patch 2
    # level 2: children of qbox 1 of patch 1 and qbox 2 of patch 2

    # children of qbox 2 of patch 2 on level 1
    children_lvl2 = Geometry.get_child_qbox_ids(1, 2, qbg_2D_mp, 2)
    @test length(children_lvl2) == 4
    @test children_lvl2 == [3,4,9,10]

    # hier_id = 1 is now the first active element on level 1
    qid, lvl, pid = Geometry.get_qbox_id_hier(1, qbg_2D_mp)
    @test lvl == 1
    @test qid == 2
    @test pid == 1

    qid, lvl, pid = Geometry.get_qbox_id_hier(6, qbg_2D_mp)
    @test lvl == 1
    @test qid == 1
    @test pid == 2

    qid, lvl, pid = Geometry.get_qbox_id_hier(17, qbg_2D_mp)
    @test lvl == 1
    @test qid == 6
    @test pid == 2

    qid, lvl, pid = Geometry.get_qbox_id_hier(34, qbg_2D_mp)
    @test lvl == 2
    @test qid == 5
    @test pid == 1

    qid, lvl, pid = Geometry.get_qbox_id_hier(55, qbg_2D_mp)
    @test lvl == 2
    @test qid == 10
    @test pid == 2

    qid = Geometry.get_qbox_id_local(1,6,2,qbg_2D_mp)
    @test qid == 3
end
#hier_g_2D = Geometry.get_hierarchical_geometry(qbg_2D)
#Plot.plot(hier_g_2D; vtk_filename="2D 4x4 square after QBox")

#hier_g_2D_mp = Geometry.get_hierarchical_geometry(qbg_2D_mp)
#Plot.plot(qbg_2D_mp; vtk_filename="test 2D multi-patch after QBox")

@testset "QBoxGeometry basic tests 2D Multi-patch from exisiting" begin
     # element 1 → qbox 1
    qid, lvl, pid = Geometry.get_qbox_id_hier(1, qbg_2D_mp_e)
    @test lvl == 1
    @test qid == 1
    @test pid == 1
    

    # element 7 (row 2, col 3) → qbox 2
    qid, lvl, pid = Geometry.get_qbox_id_hier(16, qbg_2D_mp_e)
    @test lvl == 1
    @test qid == 4
    @test pid == 1

    # element 28  → qbox 4
    qid, lvl, pid = Geometry.get_qbox_id_hier(38, qbg_2D_mp_e)
    @test lvl == 1
    @test qid == 2
    @test pid == 2

    # qbox 4 of patch 2 contains 4 elements on level 1 (not that this returns level ids)
    @test Geometry.get_qbox_element_ids(1, 2, qbg_2D_mp_e, 3) == [31,32,33,40,41,42]

    Geometry.refine_qbox!(qbg_2D_mp_e, 1, 1, 3)
    Geometry.refine_qbox!(qbg_2D_mp_e, 1, 2, 2)
    # Now active elements:
    # level 1: elements not in qbox 3 of patch 1 and qbox 2 of patch 2
    # level 2: children of qbox 3 of patch 1 and qbox 2 of patch 2

    # # children of qbox 2 of patch 2 on level 1
    children_lvl2 = Geometry.get_child_qbox_ids(1, 2, qbg_2D_mp_e, 2)
    @test children_lvl2 == [3,4,9,10,15,16]

    # hier_id = 1 is now the first active element on level 1
    qid, lvl, pid = Geometry.get_qbox_id_hier(1, qbg_2D_mp_e)
    @test lvl == 1
    @test qid == 1
    @test pid == 1

    qid, lvl, pid = Geometry.get_qbox_id_hier(13, qbg_2D_mp_e)
    @test lvl == 1
    @test qid == 4
    @test pid == 1

    qid, lvl, pid = Geometry.get_qbox_id_hier(24, qbg_2D_mp_e)
    @test lvl == 1
    @test qid == 3
    @test pid == 2

    qid, lvl, pid = Geometry.get_qbox_id_hier(33, qbg_2D_mp_e)
    @test lvl == 2
    @test qid == 13
    @test pid == 1

    for child in Geometry.get_child_qbox_ids(1, 2, qbg_2D_mp_e, 2)
        println(child, " → ", Geometry.get_qbox_element_ids(2, 2, qbg_2D_mp_e, child))
    end

    qid, lvl, pid = Geometry.get_qbox_id_hier(70, qbg_2D_mp_e)
    @test lvl == 2
    @test qid == 4
    @test pid == 2
end

@testset "QBoxGeometry basic tests 3D Multi-patch from existing" begin
    @test Geometry.get_n_elements_patch_dim(geom_3D_mp, 1) == (2,2,2)
    @test Geometry.get_n_elements_patch_dim(geom_3D_mp, 2) == (3,1,2)
    
    qid, lvl, pid = Geometry.get_qbox_id_hier(1, qbg_3D_mp)
    @test (qid, lvl, pid) == (1, 1, 1)

    qid, lvl, pid = Geometry.get_qbox_id_hier(8, qbg_3D_mp)
    @test (qid, lvl, pid) == (8, 1, 1)

    qid, lvl, pid = Geometry.get_qbox_id_hier(9, qbg_3D_mp)
    @test (qid, lvl, pid) == (1, 1, 2)

    qid, lvl, pid = Geometry.get_qbox_id_hier(14, qbg_3D_mp)
    @test (qid, lvl, pid) == (6, 1, 2)

    Geometry.refine_qbox!(qbg_3D_mp, 1, 1, 1)
    Geometry.refine_qbox!(qbg_3D_mp, 1, 2, 3)

    active = Geometry.get_active_elements(Geometry.get_hierarchical_geometry(qbg_3D_mp))
    @show Hierarchy.get_level_ids(active)

    children = Geometry.get_child_qbox_ids(1, 2, qbg_3D_mp, 3)
    @test length(children) == 8
    @test children == [5,6,11,12,17,18,23,24]

    for (i, child) in enumerate(children)
        elems = Geometry.get_qbox_element_ids(2, 2, qbg_3D_mp, child)
        @test length(elems) == 1   #  qbox_size = (1,1,1)
    end

    qid, lvl, pid = Geometry.get_qbox_id_hier(21, qbg_3D_mp)
    @test lvl == 2
    @test pid == 2
    @test qid == 5

    qid, lvl, pid = Geometry.get_qbox_id_hier(24, qbg_3D_mp)
    @test lvl == 2
    @test pid == 2
    @test qid == 12

    qid, lvl, pid = Geometry.get_qbox_id_hier(14, qbg_3D_mp)
    @test lvl == 2
    @test pid == 1
    @test qid == 2
end
#Plot.plot(qbg_3D_mp; vtk_filename="Basic 3D multi-patch_e after QBox")
end