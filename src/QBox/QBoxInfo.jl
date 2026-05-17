using ..Points
using ..Hierarchy
using ..Geometry

############################################################################################
#                                        Structure                                         #
############################################################################################

"""
        struct QBoxInfo

    Contains information about the spline degrees and number of elements of a mesh, which is
    the basis for the PBox structure. The spline degrees are constructed in such a way that
    the degree in dimension i on level 'l-1' is smaller than or equal to the degree in
    dimension i on level l.
"""
struct QBoxInfo{manifold_dim}
    qbox_size::NTuple{manifold_dim,Int}
    n_qboxes_level1::NTuple{manifold_dim,Int}
    refinement_factors::Vector{NTuple{manifold_dim,Int}}
    geometry::AbstractGeometry
    active_info::ActiveInfo
    # function QBoxInfo(qbox_size, n_qboxes_level1)
    #     n_elements = qbox_size .* n_qboxes_level1
        

    #     # for i in 1:(length(p_per_level)-1)
    #     #     if any(p_per_level[i] .> p_per_level[i+1])
    #     #         throw(ArgumentError("Spline degrees are not valid, they must be non-decreasing per level."))
    #     # end
    #     return new(n_elements_level1, p_per_level)
    # end
end

############################################################################################
#                                         Getters                                          #
############################################################################################

function get_n_qboxes_level1(qbox_info::QBoxInfo)
    return qbox_info.n_qboxes_level1
end

function get_qbox_size(qbox_info::QBoxInfo)
    return qbox_info.qbox_size
end

function get_qbox_geometry(qbox_info::QBoxInfo)
    return qbox_info.geometry
end

function get_qbox_size(qbox_info::QBoxInfo, manifold_dim::Int)
    return get_qbox_size(qbox_info)[manifold_dim]
end

function get_qbox_active_info(qbox_info::QBoxInfo)
    return qbox_info.active_info
end

function get_refinement_factors(qbox_info::QBoxInfo)
    return qbox_info.refinement_factors
end

function get_refinement_factors(qbox_info::QBoxInfo, level::Int)
    return qbox_info.refinement_factors[level]
end

"""
get_n_qboxes_dim(qbox_info::QBoxInfo, level::Int)
Returns a tuple with the number of qboxes per dimension of a certain level
"""
function get_n_qboxes_dim(qbox_info::QBoxInfo, level::Int)
    n_qboxes = get_n_qboxes_level1(qbox_info)
    refinement = get_refinement_factors(qbox_info)
    if level == 1
        return n_qboxes
    end
    for rf in refinement[1:(level-1)]
        n_qboxes = ntuple(i -> n_qboxes[i] * rf[i], length(n_qboxes))
    end
    return n_qboxes
end

"""
get_n_qboxes_total(qbox_info::QBoxInfo, level::Int)
returns the total number of qboxes in on that level
"""
function get_n_qboxes_total(qbox_info::QBoxInfo, level::Int)
    return prod(get_n_qboxes_dim(qbox_info, level))
end

"""
get_n_elements_dim(qbox_info::QBoxInfo, level::Int)
Returns a tuple with the number of elements per dimension of a certain level
"""
function get_n_elements_dim(qbox_info::QBoxInfo, level::Int)
    return get_n_qboxes_dim(qbox_info, level) .* get_qbox_size(qbox_info)
end

"""
get_n_elements_total(qbox_info::QBoxInfo, level::Int)
returns the total number of elements in on that level
"""
function get_n_elements_total(qbox_info::QBoxInfo, level::Int)
    return prod(get_n_elements_dim(qbox_info, level))
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
function get_qbox_id_hier(hier_id::Int, qbox_info::QBoxInfo)
    # hier_id →  convert_to_level_and_level_id  →  level + level_id
    level, level_id = convert_to_level_and_level_id(get_qbox_active_info(qbox_info), hier_id)

    # level_id →  get_patch_and_local_element_id  →  patch_id + local_id
    patch_id, local_element_id = get_patch_and_local_element_id(get_qbox_geometry(qbox_info), level_id)

    qbox_id= get_qbox_id_local(local_element_id, qbox_info, level)

    return qbox_id, level, patch_id
end

function get_qbox_id_local(local_element_id::Int, qbox_info::QBoxInfo, level::Int)
    # find qbox size and n_elements of the correct level
    size_qbox = get_qbox_size(qbox_info)  # NTuple{manifold_dim,Int}
    n_elements = get_n_of_elements_of_level(qbox_info, level)

    # qbox_size + n_elements + local_id → qbox_id
    coords = CartesianIndices(n_elements)[local_element_id]
    qbox_coords = ntuple(i -> (coords[i] - 1) ÷ size_qbox[i] + 1, length(size_qbox))
    n_qboxes = get_n_qboxes_dim(qbox_info, level) 
    qbox_id = LinearIndices(n_qboxes)[qbox_coords]
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
    n_elements = get_n_of_elements_of_level(qbox_info, level)
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
