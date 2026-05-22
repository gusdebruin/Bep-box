############################################################################################
#                                        Structure                                         #
############################################################################################
struct QBoxGeometry{manifold_dim, image_dim, num_patches, HG} <:
       AbstractGeometry{manifold_dim, image_dim, num_patches}
    hier_geom::HG
    qbox_size::NTuple{manifold_dim,Int}

    function QBoxGeometry(
        hier_geom::HG, qbox_size::NTuple{manifold_dim,Int}
    ) where {
        manifold_dim,
        image_dim,
        num_patches,
        HG <: HierarchicalGeometry{manifold_dim, image_dim, num_patches},
    }
        total = Hierarchy.get_num_objects(get_active_elements(hier_geom))
        if total % prod(qbox_size) != 0
            throw(
                ArgumentError(
                    "Number of 'active_elements' must be divisible by the 'qbox_size'." 
                ),
            )
        end
        return new{manifold_dim, image_dim, num_patches, HG}(hier_geom, qbox_size)
    end
end

function refine_geometry(geom::AbstractGeometry{manifold_dim, image_dim, num_patches}, 
    qbox_size::NTuple{manifold_dim,Int}) where {manifold_dim, image_dim, num_patches}
    patches = get_num_patches(geom)
    refined_patches = Vector{AbstractGeometry}(undef, patches)
    i=1
    n_elements = 0
    for patch in 1:patches
        patch_geom = get_parametric_geometry(geom, patch)
        refined_geom = subdivide_geometry(patch_geom, qbox_size)
        n_elements +=get_lin_num_elements(refined_geom)
        refined_patches[i]=refined_geom
        i +=1
    end
    return refined_patches, n_elements
end

function QBoxGeometry(geom::CartesianGeometry{manifold_dim, image_dim, num_patches}, 
    qbox_size::NTuple{manifold_dim,Int}) where {manifold_dim, image_dim, num_patches}
    refined_patches, n_elements = refine_geometry(geom, qbox_size)
    refined_geom = MultiPatchGeometry((refined_patches...))
    active_elements = Hierarchy.ActiveInfo([collect(1:n_elements)])
    hier_geom = HierarchicalGeometry((refined_geom,), active_elements)
    return QBoxGeometry(hier_geom, qbox_size)
end

function QBoxGeometry(geom::MappedGeometry{manifold_dim, image_dim, num_patches}, 
    qbox_size::NTuple{manifold_dim,Int}) where {manifold_dim, image_dim, num_patches}
    mapping = get_mapping(geom)
    refined_patches, n_elements = refine_geometry(geom, qbox_size)
    refined_geom = MappedGeometry(MultiPatchGeometry((refined_patches...)), mapping)
    active_elements = Hierarchy.ActiveInfo([collect(1:n_elements)])
    hier_geom = HierarchicalGeometry((refined_geom,), active_elements)
    return QBoxGeometry(hier_geom, qbox_size)
end

function QBoxGeometry(geom::MaskedGeometry{manifold_dim, image_dim, num_patches}, 
    qbox_size::NTuple{manifold_dim,Int}) where {manifold_dim, image_dim, num_patches}
    mask = get_evaluation_mask(geom)
    refined_patches, n_elements = refine_geometry(geom, qbox_size)
    refined_geom = MaskedGeometry(MultiPatchGeometry((refined_patches...)), mask)
    active_elements = Hierarchy.ActiveInfo([collect(1:n_elements)])
    hier_geom = HierarchicalGeometry((refined_geom,), active_elements)
    return QBoxGeometry(hier_geom, qbox_size)
end

############################################################################################
#                                         Getters                                          #
############################################################################################

function get_hierarchical_geometry(qbox_geometry::QBoxGeometry)
    return qbox_geometry.hier_geom
end

function get_qbox_size(qbox_geometry::QBoxGeometry)
    return qbox_geometry.qbox_size
end

# TODO: Check if I want to use this function; if so, replace the lines of code in other functions.
function get_patch_elements_dim(qbox_geometry::QBoxGeometry, level::Int, patch_id::Int)
    hier_geom = get_hierarchical_geometry(qbox_geometry)
    geom = hier_geom.geometries[level]
    patch_geom = get_parametric_geometry(geom, patch_id)
    return get_cart_num_elements(patch_geom)
end

function get_qbox_id_hier(hier_id::Int, qbox_geometry::QBoxGeometry)
    hier_geom = get_hierarchical_geometry(qbox_geometry)
    # hier_id →  convert_to_level_and_level_id (from ActiveInfo --> module Hierarchy)  →  level + level_id
    level, level_id = Hierarchy.convert_to_level_and_level_id(get_active_elements(hier_geom), hier_id)

    # level_id →  get_patch_and_local_element_id (from module Geometry) →  patch_id + local_id
    geom_level = get_level_geometry(hier_geom, level)
    patch_id, local_element_id = get_patch_and_local_element_id(geom_level, level_id)

    qbox_id= get_qbox_id_local(level, local_element_id, patch_id, qbox_geometry)

    return qbox_id, level, patch_id
end

function get_qbox_id_local(level::Int, local_element_id::Int, patch_id::Int, qbox_geometry::QBoxGeometry)
    # find qbox size and n_elements of the correct level
    size_qbox = get_qbox_size(qbox_geometry)
    hier_geom = get_hierarchical_geometry(qbox_geometry)
    geom_level = get_level_geometry(hier_geom, level)
    patch_geom = get_parametric_geometry(geom_level, patch_id)
    n_elements_dim = get_cart_num_elements(patch_geom)

    # qbox_size + n_elements + local_id → qbox_id
    coords = Points.CartesianIndices(n_elements_dim)[local_element_id]
    qbox_coords = ntuple(i -> (coords[i] - 1) ÷ size_qbox[i] + 1, length(size_qbox))
    n_qboxes_dim = ntuple(i -> n_elements_dim[i] ÷ size_qbox[i], length(size_qbox))
    qbox_id = Points.LinearIndices(n_qboxes_dim)[qbox_coords]
    return qbox_id
end

function get_qbox_element_ids(level::Int, patch_id::Int, qbox_geometry::QBoxGeometry, qbox_id::Int)
    size_qbox = get_qbox_size(qbox_geometry)
    hier_geom = get_hierarchical_geometry(qbox_geometry)
    geom_level = get_level_geometry(hier_geom, level)
    patch_geom = get_parametric_geometry(geom_level, patch_id)
    n_elements_dim = get_cart_num_elements(patch_geom)
    n_qboxes_dim = ntuple(i -> n_elements_dim[i] ÷ size_qbox[i], length(size_qbox))
    qbox_coords = Points.CartesianIndices(n_qboxes_dim)[qbox_id]
    
    n = prod(size_qbox)
    ids = Vector{Int}(undef, n)
    k=1
    element_ranges = ntuple(i -> (qbox_coords[i]-1)*size_qbox[i]+1 : qbox_coords[i]*size_qbox[i], length(size_qbox))
    for element_coords in Iterators.product(element_ranges...)
        local_element_id = Points.LinearIndices(n_elements_dim)[element_coords...]
        level_element_id = get_global_element_id(geom_level, patch_id, local_element_id)
        ids[k] = level_element_id
        k += 1
    end
    return ids
end

#TODO: Look at check, ask supervisors if it is good enough
function get_child_qbox_ids(level::Int, patch_id::Int, qbox_geometry::QBoxGeometry, qbox_id::Int)
    size_qbox = get_qbox_size(qbox_geometry)
    hier_geom = get_hierarchical_geometry(qbox_geometry)

    # check if level isn't the last level
    last_level = get_num_levels(hier_geom)
    if level == last_level
        throw(
                ArgumentError(
                    "Level of QBox is last level, create new geometry for new level"
                ),
            )
    end

    geom_parent = get_level_geometry(hier_geom, level)
    patch_geom_parent = get_parametric_geometry(geom_parent, patch_id)
    n_elements_dim_parent = get_cart_num_elements(patch_geom_parent)

    geom_child = get_level_geometry(hier_geom, level+1)
    patch_geom_child = get_parametric_geometry(geom_child, patch_id)
    n_elements_dim_child = get_cart_num_elements(patch_geom_child)

    rf = ntuple(i -> n_elements_dim_child[i] ÷ n_elements_dim_parent[i], length(n_elements_dim_parent))

    n_qboxes_dim_parent = ntuple(i -> n_elements_dim_parent[i] ÷ size_qbox[i], length(size_qbox))
    n_qboxes_dim_child = ntuple(i -> n_elements_dim_child[i] ÷ size_qbox[i], length(size_qbox))

    parent_coords = Points.CartesianIndices(n_qboxes_dim_parent)[qbox_id]
    child_ranges = ntuple(i -> ((parent_coords[i]-1)*rf[i] + 1) : (parent_coords[i]*rf[i]), length(rf))

    children = Int[]
    for child_coords in Iterators.product(child_ranges...)
        child_id = LinearIndices(n_qboxes_dim_child)[child_coords...]
        push!(children, child_id)
    end
    return children
end

############################################################################################
#                                       Refinement                                         #
############################################################################################

"""
QBox refinement works per patch; QBoxes are patch-local.
"""
function refine_qbox!(qbox_geometry::QBoxGeometry, level::Int, patch_id::Int, qbox_id::Int)
    hier_geom = get_hierarchical_geometry(qbox_geometry)
    active = get_active_elements(hier_geom)
    
    remove = get_qbox_element_ids(level, patch_id, qbox_geometry, qbox_id)
    
    children = get_child_qbox_ids(level, patch_id, qbox_geometry, qbox_id)
    add = Int[]  
    for child in children
        append!(add, get_qbox_element_ids(level+1, patch_id, qbox_geometry, child))
    end
    update!(active, level, remove, add)
end
