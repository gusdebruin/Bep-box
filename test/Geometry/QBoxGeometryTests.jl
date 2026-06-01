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
# QBox size
qbox_size_1D = (2,)

# Number of subdivisions
num_subdivisions_1D = (2,)
# # Level 2: refined geometry
# geom_lvl2_1D, n_elements_lvl2_1D = Geometry.refine_geometry(geom_lvl1_1D, qbox_size_1D)

# # ActiveInfo: all elements of level 1 are active
# active_1D = Hierarchy.ActiveInfo([
#     collect(1:n_elements_lvl1_1D),  
#     Int[]           
# ])
# #create the hierarchical geometry and a QBox geometry
# hier_1D = Geometry.HierarchicalGeometry((geom_lvl1_1D, geom_lvl2_1D), active_1D)
qbg_1D = Geometry.QBoxGeometry(geom_lvl1_1D, qbox_size_1D, num_subdivisions_1D)

#2D QBox three-level hierarchy
geom_lvl1_2D = Geometry.CartesianGeometry((
    0.0:2.0:4.0,
    0.0:2.0:4.0
))
n1_2D = Geometry.get_num_elements(geom_lvl1_2D)   # = 16
qbox_size_2D = (2,2)
num_subdivisions_2D = (2,2)



# geom_lvl2_2D, n2_2D = Geometry.refine_geometry(geom_lvl1_2D, qbox_size_2D)
# geom_lvl3_2D, n3_2D = Geometry.refine_geometry(geom_lvl2_2D, qbox_size_2D)

# # ActiveInfo: all elements of level 1 are active
# active_2D = Hierarchy.ActiveInfo([collect(1:n1_2D), Int[], Int[]])

# hier_2D = Geometry.HierarchicalGeometry(
#     (geom_lvl1_2D, geom_lvl2_2D, geom_lvl3_2D),
#     active_2D
# )
qbg_2D = Geometry.QBoxGeometry(geom_lvl1_2D, qbox_size_2D, num_subdivisions_2D)

# 2D multi patch geometry
geom_lvl1_2D_mp = Geometry.CartesianGeometry(((LinRange(0.5, 2.5, 3), LinRange(-0.75, 0.75, 2)),(LinRange(2.5, 5.5, 4), LinRange(-0.75, 0.75, 3))))

n1_2D_mp = Geometry.get_num_elements(geom_lvl1_2D_mp) 
qbox_size_2D_mp = (2,2)
num_subdivisions_2D_mp = (2,2)

# geom_lvl2_2D_mp, n2_2D_mp = Geometry.refine_geometry(geom_lvl1_2D_mp, qbox_size_2D_mp)
# #geom_lvl3_2D_mp, n3_2D_mp = Geometry.refine_geometry(geom_lvl2_2D_mp, qbox_size_2D_mp)

# # ActiveInfo: all elements of level 1 are active
# active_2D_mp = Hierarchy.ActiveInfo([collect(1:n1_2D_mp), Int[]])

# hier_2D_mp = Geometry.HierarchicalGeometry(
#     (geom_lvl1_2D_mp, geom_lvl2_2D_mp),
#     active_2D_mp
# )
qbg_2D_mp = Geometry.QBoxGeometry(geom_lvl1_2D_mp, qbox_size_2D_mp, num_subdivisions_2D_mp)
############################################################################################
#                                       Basic Tests                                        #
############################################################################################

@testset "QBoxGeometry basic tests 1D" begin
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



@testset "QBoxGeometry basic tests 2D" begin
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

@testset "QBoxGeometry basic tests 2D Multi-patch" begin
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
end