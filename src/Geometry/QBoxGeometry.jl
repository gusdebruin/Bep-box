############################################################################################
#                                        Structure                                         #
############################################################################################
mutable struct QBoxGeometry{manifold_dim, image_dim, num_patches} <:
       AbstractGeometry{manifold_dim, image_dim, num_patches}
    hier_geom::HierarchicalGeometry{manifold_dim, image_dim, num_patches}
    qbox_size::NTuple{manifold_dim,Int}
    num_subdivisions::NTuple{manifold_dim,Int}

end

function QBoxGeometry_refine(geom::AbstractGeometry{manifold_dim, image_dim, num_patches}, 
    qbox_size::NTuple{manifold_dim,Int}, num_subdivisions::NTuple{manifold_dim,Int}) where {manifold_dim, image_dim, num_patches}
    refined_geom = subdivide_geometry(geom, qbox_size)
    n_elements = get_num_elements(refined_geom)
    active_elements = Hierarchy.ActiveInfo([collect(1:n_elements)])
    hier_geom = HierarchicalGeometry((refined_geom,), active_elements)
    return QBoxGeometry(hier_geom, qbox_size, num_subdivisions)
end

function QBoxGeometry_from_existing(geom::AbstractGeometry{manifold_dim, image_dim, num_patches}, 
    qbox_size::NTuple{manifold_dim,Int}, num_subdivisions::NTuple{manifold_dim,Int}) where {manifold_dim, image_dim, num_patches}
    for patch_id in 1:get_num_patches(geom)
        size =  get_n_elements_patch_dim(geom, patch_id)
        if any(size .% qbox_size .!= 0)
            throw(ArgumentError("Geometry is not divisible by qbox_size on patch $patch_id")) 
        end
    end
    n_elements = get_num_elements(geom)
    active_elements = Hierarchy.ActiveInfo([collect(1:n_elements)])
    hier_geom = HierarchicalGeometry((geom,), active_elements)
    return QBoxGeometry(hier_geom, qbox_size, num_subdivisions)
end

function QBoxGeometry_refine(
    geom::HierarchicalGeometry{manifold_dim, image_dim, num_patches},
    qbox_size::NTuple{manifold_dim,Int}, 
    num_subdivisions::NTuple{manifold_dim,Int}
)where {manifold_dim, image_dim, num_patches}
    throw(ArgumentError("Cannot construct QBoxGeometry from HierarchicalGeometry"))
end

function QBoxGeometry_refine(
    geom::QBoxGeometry{manifold_dim, image_dim, num_patches},
    qbox_size::NTuple{manifold_dim,Int}, 
    num_subdivisions::NTuple{manifold_dim,Int}
)where {manifold_dim, image_dim, num_patches}
    throw(ArgumentError("Cannot construct QBoxGeometry from QBoxGeometry"))
end

function QBoxGeometry_from_existing(
    geom::HierarchicalGeometry{manifold_dim, image_dim, num_patches},
    qbox_size::NTuple{manifold_dim,Int}, 
    num_subdivisions::NTuple{manifold_dim,Int}
)where {manifold_dim, image_dim, num_patches}
    throw(ArgumentError("Cannot construct QBoxGeometry from HierarchicalGeometry"))
end

function QBoxGeometry_from_existing(
    geom::QBoxGeometry{manifold_dim, image_dim, num_patches},
    qbox_size::NTuple{manifold_dim,Int}, 
    num_subdivisions::NTuple{manifold_dim,Int}
)where {manifold_dim, image_dim, num_patches}
    throw(ArgumentError("Cannot construct QBoxGeometry from QBoxGeometry"))
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

function get_num_subdivisions(qbox_geometry::QBoxGeometry)
    return qbox_geometry.num_subdivisions
end

function get_num_elements(geometry::QBoxGeometry)
    hier_geom = get_hierarchical_geometry(geometry)
    return get_num_elements(hier_geom) 
end

function get_num_elements(geometry::QBoxGeometry, patch_id::Int)
    hier_geom = get_hierarchical_geometry(geometry)
    return get_num_elements(hier_geom, patch_id)
end

function get_num_elements_per_patch(geometry::QBoxGeometry)
    hier_geom = get_hierarchical_geometry(geometry)
    return get_num_elements_per_patch(hier_geom) 
end

function get_geometries(geometry::QBoxGeometry)
    hier_geom = get_hierarchical_geometry(geometry)
    return get_geometries(hier_geom)
end

function get_level_geometry(geometry::QBoxGeometry, level::Int)
    hier_geom = get_hierarchical_geometry(geometry)
    return get_level_geometry(hier_geom, level)
end

function get_qbox_id_hier(hier_id::Int, qbox_geometry::QBoxGeometry)
    hier_geom = get_hierarchical_geometry(qbox_geometry)
    level, level_id = Hierarchy.convert_to_level_and_level_id(get_active_elements(hier_geom), hier_id)
    
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
    n_elements_dim = get_n_elements_patch_dim(geom_level, patch_id)

    # qbox_size + n_elements + local_id → qbox_id
    coords = Points.CartesianIndices(n_elements_dim)[local_element_id]
    qbox_coords = ntuple(i -> (coords[i] - 1) ÷ size_qbox[i] + 1, length(size_qbox))
    n_qboxes_dim = ntuple(i -> n_elements_dim[i] ÷ size_qbox[i], length(size_qbox))
    lin_ind = Points.LinearIndices(n_qboxes_dim)
    qbox_id = (length(n_qboxes_dim) == 1 ? lin_ind[qbox_coords[1]] : lin_ind[qbox_coords...])
    return qbox_id
end

function add_new_level_to_geometry!(qbox_geometry::QBoxGeometry)
    hier_geo = get_hierarchical_geometry(qbox_geometry)
    num_subdivisions = get_num_subdivisions(qbox_geometry)
    new_hier_geo = subdivide_geometry(hier_geo, num_subdivisions)
    qbox_geometry.hier_geom = new_hier_geo
    return new_hier_geo
end

function get_qbox_element_ids(level::Int, patch_id::Int, qbox_geometry::QBoxGeometry, qbox_id::Int)
    size_qbox = get_qbox_size(qbox_geometry)
    hier_geom = get_hierarchical_geometry(qbox_geometry)
    geom_level = get_level_geometry(hier_geom, level)
    n_elements_dim = get_n_elements_patch_dim(geom_level, patch_id)
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

function get_child_qbox_ids(level::Int, patch_id::Int, qbox_geometry::QBoxGeometry, qbox_id::Int)
    size_qbox = get_qbox_size(qbox_geometry)
    hier_geom = get_hierarchical_geometry(qbox_geometry)

    geom_parent = get_level_geometry(hier_geom, level)
    n_elements_dim_parent = get_n_elements_patch_dim(geom_parent, patch_id)

    geom_child = get_level_geometry(hier_geom, level+1)
    n_elements_dim_child = get_n_elements_patch_dim(geom_child, patch_id)

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
    last_level = get_num_levels(hier_geom)
    if level == last_level
        hier_geom = add_new_level_to_geometry!(qbox_geometry)
    end
    active = get_active_elements(hier_geom)
    
    remove = get_qbox_element_ids(level, patch_id, qbox_geometry, qbox_id)
    
    children = get_child_qbox_ids(level, patch_id, qbox_geometry, qbox_id)
    add = Int[]  
    for child in children
        append!(add, get_qbox_element_ids(level+1, patch_id, qbox_geometry, child))
    end
    Hierarchy.update!(active, level, remove, add)
    return remove
end

function refine_qboxgeom_avg!(qbox_geometry::QBoxGeometry, errors::AbstractVector{<:Real}, dorfler::Real)
    qbox_errors = Dict{Tuple{Int,Int,Int}, Vector{Float64}}()

    for hier_id in eachindex(errors)
        qbox_id, level, patch_id = get_qbox_id_hier(hier_id, qbox_geometry)
        push!(get!(qbox_errors, (qbox_id, level, patch_id), Float64[]), errors[hier_id])
    end

    qbox_avg = Dict{Tuple{Int,Int,Int}, Float64}()
    for (key, vals) in qbox_errors
        qbox_avg[key] = sum(vals)/length(vals)
    end

    max_val, _= findmax(qbox_avg)
    threshold = (1-dorfler)*max_val
    marked_qboxes = [key for (key, avg) in qbox_avg if avg >= threshold]
    elements_that_will_be_refined = Int[]
    for (qbox_id, level, patch_id) in marked_qboxes
        remove = refine_qbox!(qbox_geometry, level, patch_id, qbox_id)
        union!(elements_that_will_be_refined, remove)
    end
    return elements_that_will_be_refined
end


function refine_qboxgeom_max!(qbox_geometry::QBoxGeometry, errors::AbstractVector{<:Real}, dorfler::Real)
    max_val,_ = findmax(errors)
    threshold = (1-dorfler)*max_val

    marked_hier_ids = findall(e -> e ≥ threshold, errors)

    marked_qboxes = Set{Tuple{Int,Int,Int}}()
    for hier_id in marked_hier_ids
        qbox_id, level, patch_id = get_qbox_id_hier(hier_id, qbox_geometry)
        push!(marked_qboxes, (qbox_id, level, patch_id))
    end

    elements_that_will_be_refined = Int[]
    for (qbox_id, level, patch_id) in marked_qboxes
        remove = refine_qbox!(qbox_geometry, level, patch_id, qbox_id)
        union!(elements_that_will_be_refined, remove)
    end
    return elements_that_will_be_refined
end