module QBoxGeometryTests

using Mantis
using Test

############################################################################################
#                                          Setup                                           #
############################################################################################

# 1D QBox two-level hierarchy
# create a 1D geometry with 4 elements : [[0,1],[1,2],[2,3],[3,4]] for level 1
geom_lvl1_1D = Geometry.CartesianGeometry((0.0:1.0:4.0,))
n_elements_lvl1_1D = Geometry.get_num_elements(geom_lvl1_1D)
# QBox size
qbox_size_1D = (2,)
# Level 2: refined geometry
geom_lvl2_1D, n_elements_lvl2_1D = Geometry.refine_geometry(geom_lvl1_1D, qbox_size_1D)

# ActiveInfo: all elements of level 1 are active
active_1D = Hierarchy.ActiveInfo([
    collect(1:n_elements_lvl1_1D),  
    Int[]           
])
#create the hierarchical geometry and a QBox geometry
hier_1D = Geometry.HierarchicalGeometry((geom_lvl1_1D, geom_lvl2_1D), active_1D)
qbg_1D = Geometry.QBoxGeometry(hier_1D, qbox_size_1D)

#2D QBox three-level hierarchy
geom_lvl1_2D = Geometry.CartesianGeometry((
    0.0:1.0:4.0,
    0.0:1.0:4.0
))
n1_2D = Geometry.get_num_elements(geom_lvl1_2D)   # = 16
qbox_size_2D = (2,2)

geom_lvl2_2D, n2_2D = Geometry.refine_geometry(geom_lvl1_2D, qbox_size_2D)
geom_lvl3_2D, n3_2D = Geometry.refine_geometry(geom_lvl2_2D, qbox_size_2D)

# ActiveInfo: all elements of level 1 are active
active_2D = Hierarchy.ActiveInfo([collect(1:n1_2D), Int[], Int[]])

hier_2D = Geometry.HierarchicalGeometry(
    (geom_lvl1_2D, geom_lvl2_2D, geom_lvl3_2D),
    active_2D
)
qbg_2D = Geometry.QBoxGeometry(hier_2D, qbox_size_2D)


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

    # children of qbox 1 (level 1 → level 2)
    children = Geometry.get_child_qbox_ids(1, 1, qbg_1D, 1)
    @test length(children) == 2

    # refine qbox 1
    Geometry.refine_qbox!(qbg_1D, 1, 1, 1)

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
    @test n2_2D == 64 
    @test n3_2D == 256

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

    # children of qbox 1 on level 2
    children_lvl2 = Geometry.get_child_qbox_ids(1, 1, qbg_2D, 1)
    @test length(children_lvl2) == 4

    Geometry.refine_qbox!(qbg_2D, 1, 1, 1)
    # Now active elements:
    # level 1: elements not in qbox 1
    # level 2: children of qbox 1

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

end