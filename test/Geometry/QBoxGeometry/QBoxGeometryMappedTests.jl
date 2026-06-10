module QBoxGeometryMappedTests
using Mantis
using Test

############################################################################################
#                                          Setup                                           #
############################################################################################

# Test MappedCartesianGeometry ------------------------------------------------

# Mappings to create the deformed geometries. The mappings are defined with reference
# to the unit square [0,1]x[0,1] as parametric domain.
function mapping_patch_1_slant(x::AbstractVector{Float64}, slant_factor=0.25)
    return [x[1] + slant_factor * x[1] * x[2], x[2]]
end
function dmapping_patch_1_slant(x::AbstractVector{Float64}, slant_factor=0.25)
    return [
        [1.0 + slant_factor * x[2] slant_factor * x[1]]
        [0.0 1.0]
    ]
end
function ddmapping_patch_1_slant(x::AbstractVector{Float64}, slant_factor=0.25)
    return (
        [
            [0.0 slant_factor]
            [slant_factor 0.0]
        ],
        [
            [0.0 0.0]
            [0.0 0.0]
        ],
    )
end
mapping_patch_1_slanted = Geometry.Mapping(
    (2, 2), mapping_patch_1_slant, dmapping_patch_1_slant, ddmapping_patch_1_slant
)
function mapping_patch_2_slant(x::AbstractVector{Float64}, slant_factor=0.25)
    return [x[1] + 1.0 + slant_factor * (1.0 - x[1]) * x[2], x[2]]
end
function dmapping_patch_2_slant(x::AbstractVector{Float64}, slant_factor=0.25)
    return [
        [1.0 - slant_factor * x[2] slant_factor * (1.0 - x[1])]
        [0.0 1.0]
    ]
end
function ddmapping_patch_2_slant(x::AbstractVector{Float64}, slant_factor=0.25)
    return (
        [
            [0.0 -slant_factor]
            [-slant_factor 0.0]
        ],
        [
            [0.0 0.0]
            [0.0 0.0]
        ],
    )
end
mapping_patch_2_slanted = Geometry.Mapping(
    (2, 2), mapping_patch_2_slant, dmapping_patch_2_slant, ddmapping_patch_1_slant
)
num_elements_per_dim_per_patch = ((4, 4), (5, 6))
geom_cart_patch_1 = Geometry.CartesianGeometry((
    0.0:(1.0 / num_elements_per_dim_per_patch[1][1]):1.0,
    0.0:(1.0 / num_elements_per_dim_per_patch[1][2]):1.0,
))
geom_cart_patch_2 = Geometry.CartesianGeometry((
    0.0:(1.0 / num_elements_per_dim_per_patch[2][1]):1.0,
    0.0:(1.0 / num_elements_per_dim_per_patch[2][2]):1.0,
))

# Reduction test
function mapping_I(x::AbstractVector{Float64}, slant_factor=0.25)
    return [x[1], x[2]]
end
function dmapping_I(x::AbstractVector{Float64}, slant_factor=0.25)
    return [
        [0.0 0.0]
        [0.0 0.0]
    ]
end
mapping_I_obj = Geometry.Mapping((1, 1), mapping_I, dmapping_I)
geometry1 = Geometry.MappedGeometry(
    Geometry.CartesianGeometry((LinRange(0.0, 1.0, 2),)), mapping_I_obj
)

# Mapped, explicit geometry and mapping per patch.
geom_slanted_2patch = Geometry.MappedGeometry(
    (geom_cart_patch_1, geom_cart_patch_2),
    (mapping_patch_1_slanted, mapping_patch_2_slanted),
)

# Mapped, one parametric domain with multiple mappings.
geom_slanted_2patch_oneref = Geometry.MappedGeometry(
    Geometry.CartesianGeometry(((LinRange(0.0, 1.0, 5), LinRange(0.0, 1.0, 7)),)),
    (mapping_patch_1_slanted, mapping_patch_2_slanted),
)

geom_slanted_2patch_onemap = Geometry.MappedGeometry(
    (
        Geometry.CartesianGeometry(((LinRange(0.0, 0.25, 3), LinRange(0.0, 1.0, 7)),)),
        Geometry.CartesianGeometry(((LinRange(0.25, 0.5, 4), LinRange(0.0, 1.0, 7)),)),
        Geometry.CartesianGeometry(((LinRange(0.5, 0.75, 5), LinRange(0.0, 1.0, 7)),)),
        Geometry.CartesianGeometry(((LinRange(0.75, 1.0, 6), LinRange(0.0, 1.0, 7)),)),
    ),
    mapping_patch_1_slanted,
)


qbox_size_slanted_2patch_onemap = (1,3)
n_sub_slanted_2patch_onemap=(2,1)
qbg_slanted_2patch_onemap = Geometry.QBoxGeometry_refine(geom_slanted_2patch_onemap, qbox_size_slanted_2patch_onemap, n_sub_slanted_2patch_onemap)
#Plot.plot(qbg_slanted_2patch_onemap; vtk_filename="QBOX Mapped, multiple patches with one map")

qbox_size_slanted_2patch_oneref = (2,3)
n_sub_slanted_2patch_oneref =(2,2)
qbg_slanted_2patch_oneref = Geometry.QBoxGeometry_from_existing(geom_slanted_2patch_oneref, qbox_size_slanted_2patch_oneref, n_sub_slanted_2patch_oneref)
#Plot.plot(qbg_slanted_2patch_oneref; vtk_filename="Before QBOX refinement Mapped, multiple patches with one map")

# Mapped, one parametric domain with one mapping. --> one patch
geom_slanted_1patch = Geometry.MappedGeometry(
    Geometry.CartesianGeometry(((LinRange(0.0, 1.0, 5), LinRange(0.0, 1.0, 7)),)),
    (mapping_patch_1_slanted),
)
qbox_size_slanted_1patch = (2,3)
n_sub_slanted_1patch =(2,2)
qbg_slanted_1patch = Geometry.QBoxGeometry_from_existing(geom_slanted_1patch, qbox_size_slanted_1patch, n_sub_slanted_1patch)
#Plot.plot(qbg_slanted_1patch; vtk_filename="Before QBOX refinement Mapped, one patch")


############################################################################################
#                                             Tests                                        #
############################################################################################

@testset "QBoxGeometry Mapped geom, multiple patches with one map" begin
     # element 1 → qbox 1
    qid, lvl, pid = Geometry.get_qbox_id_hier(1, qbg_slanted_2patch_onemap)
    @test lvl == 1
    @test qid == 1
    @test pid == 1
    

    # element 61 (row 5, col 5) → qbox 19
    qid, lvl, pid = Geometry.get_qbox_id_hier(51, qbg_slanted_2patch_onemap)
    @test lvl == 1
    @test qid == 6
    @test pid == 2

    # element 28  → qbox 4
    qid, lvl, pid = Geometry.get_qbox_id_hier(197, qbg_slanted_2patch_onemap)
    @test lvl == 1
    @test qid == 15
    @test pid == 4

    # qbox 4 of patch 2 contains 4 elements on level 1 (not that this returns level ids)
    @test Geometry.get_qbox_element_ids(1, 3, qbg_slanted_2patch_onemap, 7) == [105,109,113]

    Geometry.refine_qbox!(qbg_slanted_2patch_onemap, 1, 1, 1)
    Geometry.refine_qbox!(qbg_slanted_2patch_onemap, 1, 1, 4)
    Geometry.refine_qbox!(qbg_slanted_2patch_onemap, 1, 4, 29)
    Geometry.refine_qbox!(qbg_slanted_2patch_onemap, 1, 4, 30)
    
    Plot.plot(qbg_slanted_2patch_onemap; vtk_filename="After QBOX refinement Mapped, multiple patches with one map")
    # children of qbox 2 of patch 2 on level 1
    children_lvl2 = Geometry.get_child_qbox_ids(1, 1, qbg_slanted_2patch_onemap, 4)
    @test children_lvl2 == [7,8]

    # hier_id = 1 is now the first active element on level 1
    qid, lvl, pid = Geometry.get_qbox_id_hier(1, qbg_slanted_2patch_onemap)
    @test lvl == 1
    @test qid == 2
    @test pid == 1

    qid, lvl, pid = Geometry.get_qbox_id_hier(31, qbg_slanted_2patch_onemap)
    @test lvl == 1
    @test qid == 1
    @test pid == 2

    qid, lvl, pid = Geometry.get_qbox_id_hier(240, qbg_slanted_2patch_onemap)
    @test lvl == 1
    @test qid == 28
    @test pid == 4

    qid, lvl, pid = Geometry.get_qbox_id_hier(241, qbg_slanted_2patch_onemap)
    @test lvl == 2
    @test qid == 1
    @test pid == 1

    qid, lvl, pid = Geometry.get_qbox_id_hier(254, qbg_slanted_2patch_onemap)
    @test lvl == 2
    @test qid == 57
    @test pid == 4

    @test Geometry.get_qbox_element_ids(2, 4, qbg_slanted_2patch_onemap, 60) == [484, 494, 504]

end

@testset "QBoxGeometry Mapped geom, one parametric domain with multiple mappings" begin
         # element 1 → qbox 1
    qid, lvl, pid = Geometry.get_qbox_id_hier(1, qbg_slanted_2patch_oneref)
    @test lvl == 1
    @test qid == 1
    @test pid == 1
    

    # element 61 (row 5, col 5) → qbox 19
    qid, lvl, pid = Geometry.get_qbox_id_hier(16, qbg_slanted_2patch_oneref)
    @test lvl == 1
    @test qid == 4
    @test pid == 1

    # element 28  → qbox 4
    qid, lvl, pid = Geometry.get_qbox_id_hier(27, qbg_slanted_2patch_oneref)
    @test lvl == 1
    @test qid == 2
    @test pid == 2

    # qbox 4 of patch 2 contains 4 elements on level 1 (not that this returns level ids)
    @test Geometry.get_qbox_element_ids(1, 2, qbg_slanted_2patch_oneref, 4) == [39,40,43,44,47,48]

    Geometry.refine_qbox!(qbg_slanted_2patch_oneref, 1, 1, 2)
    
    #Plot.plot(qbg_slanted_2patch_oneref; vtk_filename="After QBOX refinement Mapped, one parametric domain with multiple mappings")
    # children of qbox 2 of patch 2 on level 1
    children_lvl2 = Geometry.get_child_qbox_ids(1, 1, qbg_slanted_2patch_oneref, 2)
    @test children_lvl2 == [3,4,7,8]

    # hier_id = 1 is now the first active element on level 1
    qid, lvl, pid = Geometry.get_qbox_id_hier(3, qbg_slanted_2patch_oneref)
    @test lvl == 1
    @test qid == 1
    @test pid == 1

    qid, lvl, pid = Geometry.get_qbox_id_hier(42, qbg_slanted_2patch_oneref)
    @test lvl == 1
    @test qid == 4
    @test pid == 2

    qid, lvl, pid = Geometry.get_qbox_id_hier(43, qbg_slanted_2patch_oneref)
    @test lvl == 2
    @test qid == 3
    @test pid == 1

    qid, lvl, pid = Geometry.get_qbox_id_hier(60, qbg_slanted_2patch_oneref)
    @test lvl == 2
    @test qid == 7
    @test pid == 1

    @test Geometry.get_qbox_element_ids(2, 1, qbg_slanted_2patch_oneref, 3) == [5,6,13,14,21,22]
end

@testset "QBoxGeometry Mapped geom, one patch" begin
         # element 1 → qbox 1
    qid, lvl, pid = Geometry.get_qbox_id_hier(1, qbg_slanted_1patch)
    @test lvl == 1
    @test qid == 1
    @test pid == 1
    
    qid, lvl, pid = Geometry.get_qbox_id_hier(16, qbg_slanted_1patch)
    @test lvl == 1
    @test qid == 4
    @test pid == 1

    @test Geometry.get_qbox_element_ids(1, 1, qbg_slanted_1patch, 1) == [1,2,5,6,9,10]

    Geometry.refine_qbox!(qbg_slanted_1patch, 1, 1, 2)
    
    # Plot.plot(qbg_slanted_1patch; vtk_filename="After QBOX refinement Mapped, one patch")
    children_lvl2 = Geometry.get_child_qbox_ids(1, 1, qbg_slanted_1patch, 2)
    @test children_lvl2 == [3,4,7,8]

    qid, lvl, pid = Geometry.get_qbox_id_hier(3, qbg_slanted_1patch)
    @test lvl == 1
    @test qid == 1
    @test pid == 1

    qid, lvl, pid = Geometry.get_qbox_id_hier(18, qbg_slanted_1patch)
    @test lvl == 1
    @test qid == 4
    @test pid == 1

    qid, lvl, pid = Geometry.get_qbox_id_hier(31, qbg_slanted_1patch)
    @test lvl == 2
    @test qid == 7
    @test pid == 1

    @test Geometry.get_qbox_element_ids(2, 1, qbg_slanted_1patch, 3) == [5,6,13,14,21,22]
end

end
