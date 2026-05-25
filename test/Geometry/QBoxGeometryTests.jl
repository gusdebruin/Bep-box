module QBoxGeometryTests

using Mantis
using Test

############################################################################################
#                                          Setup                                           #
############################################################################################

# # 1D, 1 patch, 4 elementen
# geom_1 = Geometry.CartesianGeometry((0.0:1.0:4.0,))
# qbox_size_1 = (2,)
# qbg_1 = Geometry.QBoxGeometry(geom, qbox_size)

# # 2D, 1 patch, 4×4 elementen
# geom_2 = Geometry.CartesianGeometry((0.0:1.0:4.0, 0.0:1.0:4.0))
# qbox_size_2 = (2,2)
# qbg_2 = Geometry.QBoxGeometry(geom, qbox_size)

# 1D, 1 patch, 4 elementen
geom = Geometry.CartesianGeometry((0.0:1.0:4.0,))
qbg = Geometry.QBoxGeometry(geom, (2,))
# element 1 → qbox 1
qid, lvl, pid = Geometry.get_qbox_id_hier(1, qbg)
# refine qbox 1
Geometry.refine_qbox!(qbg, 1, 1, 1)
# children of qbox 1
children = Geometry.get_child_qbox_ids(1, 1, qbg, 1)


@testset "QBoxGeometry basic tests" begin
    @test Geometry.get_patch_elements_dim(qbg, 1, 1) == (4,)
    @test qid == 1

    # qbox 1 contains elements 1 and 2
    @test Geometry.get_qbox_element_ids(1, 1, qbg, 1) == [1,2]
    @test length(children) == 2
end

end