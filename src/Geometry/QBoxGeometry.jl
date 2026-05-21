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

function QBoxGeometry(geom::AbstractGeometry{manifold_dim, image_dim, num_patches}, 
    qbox_size::NTuple{manifold_dim,Int}) where {manifold_dim, image_dim, num_patches}
    refined_geom = subdivide_geometry(geom, qbox_size) #from which module is this??
    n_elements = get_lin_num_elements(refined_geom) #from which module is this??
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

############################################################################################
#                                    Abstract Methods                                      #
############################################################################################

"""
    get_pbox_id(hier_id::Int, pbox_info::PBoxInfo, active_info::ActiveInfo,
        geometry::AbstractGeometry)

Returns the id of the pbox a element 'hier_id' is in. It uses these steps:
1. hier_id →  convert_to_level_and_level_id  →  level + level_id
2. level_id →  get_patch_and_local_element_id  →  patch_id + local_id
3. find p and n_elements of the correct level
4. p + n_elements + local_id → pbox_id

# Arguments
- `geometry::AbstractGeometry`: The multi-patch geometry.
- `hier_id::Int`: The hierarchical element ID.
- `pbox_info::PBoxInfo`: The info of the pbox (contains spline degrees and number of 
    elements of the first level)
- `active_info::ActiveInfo`: The active elements per level.

# Returns
- `pbox_id::Int`: The pbox ID
- `level::Int`: The level
- `patch_id::Int`: The patch ID
"""
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

"""
    get_pbox_element_ids(geometry::AbstractGeometry, pbox_id::Int, pbox_info::PBoxInfo,
        level::Int, patch_id::Int)

Returns the ids of the elements of a pbox 'pbox_id'.

# Arguments
- `geometry::AbstractGeometry`: The multi-patch geometry.
- `pbox_id::Int`: The pbox ID
- `pbox_info::PBoxInfo`: The info of the pbox (contains spline degrees and number of 
    elements of the first level)
- `level::Int`: The level
- `patch_id::Int`: The patch ID.

# Returns
- `ids::Vector{Int}`: Returns a Vector{Int} with the level-local element IDs
"""
function get_qbox_element_ids(level::Int, patch_id::Int, qbox_geometry::QBoxGeometry, qbox_id::Int)
    size_qbox = get_qbox_size(qbox_geometry)
    hier_geom = get_hierarchical_geometry(qbox_geometry)
    geom_level = get_level_geometry(hier_geom, level)
    patch_geom = get_parametric_geometry(geom_level, patch_id)
    n_elements_dim = get_cart_num_elements(patch_geom)
    n_qboxes_dim = ntuple(i -> n_elements_dim[i] ÷ size_qbox[i], length(size_qbox))
    qbox_coords = Points.CartesianIndices(n_qboxes_dim)[qbox_id]
    
    # hier nog naar kijken (is misschien efficienter:)
    # n = prod(size_qbox)
    # ids = Vector{Int}(undef, n)
    ids = Int[]
    element_ranges = ntuple(i -> (qbox_coords[i]-1)*size_qbox[i]+1 : qbox_coords[i]*size_qbox[i], length(size_qbox))
    for element_coords in Iterators.product(element_ranges...)
        local_element_id = Points.LinearIndices(n_elements_dim)[element_coords...]
        level_element_id = get_global_element_id(geom_level, patch_id, local_element_id)
        push!(ids, level_element_id)
    end
    return ids
end

"""
	 child_pbox_ids(pbox_info::PBoxInfo, level::Int, pbox_id::Int)
Returns the ids of the 'children' of `pbox_id`, the children of a pbox 'pbox_id' are the 
pboxes in which 'pbox_id' has been refined to in the next level.

# Arguments
- `pbox_info::PBoxInfo`: The info of the pbox (contains spline degrees and number of 
    elements of the first level)
- `level::Int`: The level of the pbox.
- `pbox_id::Int`: The pbox ID of which we want to know the children. 

# Returns
- `children::Vector{Int}`: Returns a Vector{Int} with IDs of the children

"""
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
	refine_pbox!(geometry::AbstractGeometry, active_info::ActiveInfo, pbox_info::PBoxInfo, 
        level::Int, patch_id::Int, pbox_id::Int)

It finds the elements of the pbox 'pbox_id' that needs to be refined. And it finds all the 
elements of the children of the pbox. Then using the function update! the elements of 
the children will be set to active, and the elements of the pbox no longer.

# Arguments
- `geometry::AbstractGeometry`: The multi-patch geometry.
- `active_info::ActiveInfo`: The active elements per level.
- `pbox_info::PBoxInfo`: The info of the pbox (contains spline degrees and number of 
    elements of the first level)
- `level::Int`: The level of the pbox we want to refine.
- `patch_id::Int`: The patch ID in which the pbox is that we want to refine.
- `pbox_id::Int`: The pbox ID of the pbox that we want to refine.

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
