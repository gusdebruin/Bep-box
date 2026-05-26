module QBoxGeometryTests

using Mantis
using Test

############################################################################################
#                                          Setup                                           #
############################################################################################
# create a 1D geometry with 4 elements
# Level 1: original geometry
geom_lvl1 = Geometry.CartesianGeometry((0.0:1.0:4.0,))
n1 = Geometry.get_num_elements(geom_lvl1)
# QBox size
qbox_size = (2,)


# Level 2: refined geometry
geom_lvl2, n_elements_lvl2 = Geometry.refine_geometry(geom_lvl1, qbox_size)

# ActiveInfo: all elements of level 1 are active
active = Hierarchy.ActiveInfo([
    collect(1:n1),  
    Int[]           
])

hier = Geometry.HierarchicalGeometry((geom_lvl1, geom_lvl2), active)
qbg = Geometry.QBoxGeometry(hier, qbox_size)

@testset "QBoxGeometry basic tests 1D" begin
    # element 1 → qbox 1
    qid, lvl, pid = Geometry.get_qbox_id_hier(1, qbg)
    @test qid == 1
    @test lvl == 1
    @test pid == 1

    qid, lvl, pid = Geometry.get_qbox_id_hier(3, qbg)
    @test qid == 2
    @test lvl == 1
    @test pid == 1

    # qbox 1 contains elements 1 and 2
    @test Geometry.get_qbox_element_ids(1, 1, qbg, 1) == [1,2]

    # children of qbox 1 (level 1 → level 2)
    children = Geometry.get_child_qbox_ids(1, 1, qbg, 1)
    @test length(children) == 2

    # refine qbox 1
    Geometry.refine_qbox!(qbg, 1, 1, 1)

    qid, lvl, pid = Geometry.get_qbox_id_hier(1, qbg)
    @test qid == 2
    @test lvl == 1
    @test pid == 1

    qid, lvl, pid = Geometry.get_qbox_id_hier(3, qbg)
    @test qid == 1
    @test lvl == 2
    @test pid == 1
end

end