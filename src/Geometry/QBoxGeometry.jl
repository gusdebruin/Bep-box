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
    refined_geom = subdivide_geometry(geom, qbox_size)
    n_elements = get_lin_num_elements(refined_geom)
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
    hg = get_hierarchical_geometry(qbox_geometry)
    # hier_id →  convert_to_level_and_level_id (from ActiveInfo --> module Hierarchy)  →  level + level_id
    level, level_id = Hierarchy.convert_to_level_and_level_id(get_active_elements(hg), hier_id)

    # level_id →  get_patch_and_local_element_id (from module Geometry) →  patch_id + local_id
    patch_id, local_element_id = get_patch_and_local_element_id(get_geometries(hg)[level], level_id)

    qbox_id= get_qbox_id_local(level, local_element_id, patch_id, qbox_geometry)

    return qbox_id, level, patch_id
end

function get_qbox_id_local(level::Int, local_element_id::Int, patch_id::Int, qbox_geometry::QBoxGeometry)
    # find qbox size and n_elements of the correct level
    size_qbox = get_qbox_size(qbox_geometry)  # NTuple{manifold_dim,Int}
    hg = get_hierarchical_geometry(qbox_geometry)
    geom = get_geometries(hg)[level]
    patch_geom = get_parametric_geometry(geom, patch_id)
    n_elements_dim = get_cart_num_elements(patch_geom)
    
    # qbox_size + n_elements + local_id → qbox_id
    coords = CartesianIndices(n_elements_dim)[local_element_id]
    qbox_coords = ntuple(i -> (coords[i] - 1) ÷ size_qbox[i] + 1, length(size_qbox))
    n_qboxes_dim = ntuple(i -> n_elements_dim[i] ÷ size_qbox[i], length(size_qbox))
    qbox_id = LinearIndices(n_qboxes_dim)[qbox_coords]
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
function get_qbox_element_ids(qbox_id::Int, qbox_info::QBoxInfo, level::Int, patch_id::Int)
    q = get_qbox_size(qbox_info)    
    n_elements = get_n_elements_dim(qbox_info, level)
    n_qboxes = get_n_qboxes_dim(qbox_info, level)
    qbox_coords = CartesianIndices(n_qboxes)[qbox_id]

    ids = Int[]
    element_ranges = ntuple(i -> (qbox_coords[i]-1)*q[i]+1 : qbox_coords[i]*q[i], length(q))
    geom = get_qbox_geometry(qbox_info)
    for element_coords in Iterators.product(element_ranges...)
        element_local_id = LinearIndices(n_elements)[element_coords...]
        element_level_id = get_global_element_id(geom, patch_id, element_local_id)
        push!(ids, element_level_id)
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
function child_qbox_ids(qbox_info::QBoxInfo, level::Int, qbox_id::Int)
    refinement = get_refinement_factors(qbox_info, level)
    n_qboxes_parent = get_n_qboxes_dim(qbox_info, level)
    n_qboxes_child  = get_n_qboxes_dim(qbox_info, level + 1)
    
    qbox_coords = CartesianIndices(n_qboxes_parent)[qbox_id]

    offset_ranges = ntuple(i -> 0:(refinement[i]-1), length(refinement))
    offsets = Iterators.product(offset_ranges...)

    children = Int[]
    for off in offsets
        child_coords = ntuple(i -> (qbox_coords[i]-1)*refinement[i] + off[i] + 1, length(refinement))
        child_id = LinearIndices(n_qboxes_child)[child_coords...]
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
function refine_qbox!(qbox_info::QBoxInfo, level::Int, patch_id::Int, qbox_id::Int)
    active = get_qbox_active_info(qbox_info)
    remove = get_qbox_element_ids(qbox_id, qbox_info, level, patch_id)

    children = child_qbox_ids(qbox_info, level, qbox_id)
    add = Int[]  
    for child in children
        append!(add, get_qbox_element_ids(child, qbox_info, level+1, patch_id))
    end
    update!(active, level, remove, add)
end
