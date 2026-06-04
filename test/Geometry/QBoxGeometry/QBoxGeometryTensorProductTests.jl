module QBoxGeometryTensorProductTests
using Mantis
using Test

############################################################################################
#                                          Setup                                           #
############################################################################################

# Test Square Tensor Product Geometry -----------------------------------------
# Generate a tensor product geometry by combining two lines

# Line geometries
line_1_geometry = Geometry.create_cartesian_box((0.0,), (1.0,), (10,))
line_2_geometry = Geometry.create_cartesian_box((2.0,), (1.0,), (10,))

# Tensor product geometry
tensor_prod_geometry = Geometry.TensorProductGeometry((line_1_geometry, line_2_geometry))
#Plot.plot(tensor_prod_geometry; vtk_filename="tensor_prod_geometry")
tpg_qbox_size = (2,5)
tpg_num_sub = (3,2)
qbg_tpg = Geometry.QBoxGeometry_from_existing(tensor_prod_geometry, tpg_qbox_size, tpg_num_sub)


# Test Cylinder Tensor Product Geometry ---------------------------------------
deg = 2
nθ_elements = 4
Wt = 2.0 * pi / nθ_elements
b = FunctionSpaces.GeneralizedTrigonometric(deg, Wt)
breakpoints = collect(LinRange(0.0, nθ_elements, nθ_elements + 1))
patch = Geometry.CartesianGeometry(breakpoints)
B = FunctionSpaces.BSplineSpace(patch, b, [-1, 1, 1, 1, -1])
GB = FunctionSpaces.GTBSplineSpace((B,), [1])

# control points for geometry
# radius of cylinder is 1.0
geom_coeffs_circle = [
    +1.0 -1.0
    +1.0 +1.0
    -1.0 +1.0
    -1.0 -1.0
]
cylinder_circle_geometry = FunctionSpaces.DiscreteGeometry(GB, geom_coeffs_circle)
dx_cylinder_line = 0.1
nz_elements = 10
cylinder_line_geometry = Geometry.create_cartesian_box((0.0,), (1.0,), (nz_elements,))

# Tensor product geometry
cylinder_tensor_prod_geometry = Geometry.TensorProductGeometry((
    cylinder_circle_geometry, cylinder_line_geometry
))
#Plot.plot(cylinder_tensor_prod_geometry; vtk_filename="cylinder_tensor_prod_geometry")



# Reduction test, single-patch, single element, single geometry, 1D.
cg1 = Geometry.CartesianGeometry(([-1, 1],))
tpgeometry1 = Geometry.TensorProductGeometry((cg1,))
qbox_size_1 = (4,)
num_sub_1 = (2,)
qbp_1 = Geometry.QBoxGeometry_refine(tpgeometry1, qbox_size_1, num_sub_1)


# All cartesian. 3D: 2D (2 patches) tensored with 1D (3 patches).
cg1d = Geometry.CartesianGeometry((
    (LinRange(0.0, 1.0, 2),), (LinRange(1.0, 2.0, 3),), (LinRange(2.0, 3.0, 4),)
))
cg2d = Geometry.CartesianGeometry((
    (LinRange(0.0, 1.0, 4), LinRange(0.0, 1.0, 5)), # First patch
    (LinRange(1.0, 2.0, 6), LinRange(0.0, 1.0, 7)), # Second patch
))
tpgeometry2 = Geometry.TensorProductGeometry((cg2d, cg1d))
#Plot.plot(tpgeometry2; vtk_filename="tpgeometry2")
tpg2_qbox_size = (1,2,1)
tpg2_num_sub = (2,2,2)
qbp_tpg2 = Geometry.QBoxGeometry_from_existing(tpgeometry2, tpg2_qbox_size, tpg2_num_sub)


############################################################################################
#                                             Tests                                        #
############################################################################################

@testset "QBoxGeometry TensorProduct geom, one element" begin
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

@testset "QBoxGeometry TensorProduct geom, two lines" begin
     # element 1 → qbox 1
    qid, lvl, pid = Geometry.get_qbox_id_hier(1, qbg_tpg)
    @test lvl == 1
    @test qid == 1
    @test pid == 1
    
    qid, lvl, pid = Geometry.get_qbox_id_hier(51, qbg_tpg)
    @test lvl == 1
    @test qid == 6
    @test pid == 1

    qid, lvl, pid = Geometry.get_qbox_id_hier(100, qbg_tpg)
    @test lvl == 1
    @test qid == 10
    @test pid == 1

    @test Geometry.get_qbox_element_ids(1, 1, qbg_tpg, 7) == [53,54,63,64,73,74,83,84,93,94]

    Geometry.refine_qbox!(qbg_tpg, 1, 1, 3)
    Geometry.refine_qbox!(qbg_tpg, 1, 1, 8)

    #Plot.plot(qbg_tpg; vtk_filename="After QBOX refinement Tensorproduct, two lines")
    
    children_lvl2 = Geometry.get_child_qbox_ids(1, 1, qbg_tpg, 4)
    @test children_lvl2 == [10,11,12,25,26,27]

    qid, lvl, pid = Geometry.get_qbox_id_hier(1, qbg_tpg)
    @test lvl == 1
    @test qid == 1
    @test pid == 1

    qid, lvl, pid = Geometry.get_qbox_id_hier(80, qbg_tpg)
    @test lvl == 1
    @test qid == 10
    @test pid == 1

    qid, lvl, pid = Geometry.get_qbox_id_hier(83, qbg_tpg)
    @test lvl == 2
    @test qid == 8
    @test pid == 1


    @test Geometry.get_qbox_element_ids(2, 1, qbg_tpg, 7) == [13,14,43,44,73,74,103,104,133,134]

end

@testset "QBoxGeometry tpgeometry2" begin
     # element 1 → qbox 1
    qid, lvl, pid = Geometry.get_qbox_id_hier(1, qbp_tpg2)
    @test lvl == 1
    @test qid == 1
    @test pid == 1
    
    qid, lvl, pid = Geometry.get_qbox_id_hier(6, qbp_tpg2)
    @test lvl == 1
    @test qid == 3
    @test pid == 1

    qid, lvl, pid = Geometry.get_qbox_id_hier(28, qbp_tpg2)
    @test lvl == 1
    @test qid == 6
    @test pid == 2

    @test Geometry.get_qbox_element_ids(1, 1, qbp_tpg2, 6) == [9,12]
    @test Geometry.get_qbox_element_ids(1, 3, qbp_tpg2, 6) == [51,54]

    Geometry.refine_qbox!(qbp_tpg2, 1, 1, 1)


    Plot.plot(qbp_tpg2; vtk_filename="After QBOX refinement qbp_tpg2")
    
    children_lvl2 = Geometry.get_child_qbox_ids(1, 1, qbp_tpg2, 1)
    @test children_lvl2 == [1,2,7,8,25,26,31,32]

    @test Geometry.get_qbox_element_ids(2, 1, qbp_tpg2, 2) == [2,8]

end

end