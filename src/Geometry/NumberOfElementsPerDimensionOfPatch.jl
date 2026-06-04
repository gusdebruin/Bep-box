"""
	get_n_elements_patch_dim(geom::AbstractGeometry, patch_id::Int)

Returns a tuple with per dimension the number of elements of that patch. 

# Returns
- `n_elements_dim`: The number of elements per dimension corresponding to `patch_id`.
"""
function get_n_elements_patch_dim(geom::AbstractGeometry, patch_id::Int)
    throw(MethodError(get_n_elements_patch_dim, (geom, patch_id)))
end

function get_n_elements_patch_dim(geom::CartesianGeometry, patch_id::Int)
    patch_break = get_breakpoints(geom, patch_id)
    n_elements_dim = ntuple(d -> length(patch_break[d]) - 1, length(patch_break))
    return n_elements_dim
end

function get_n_elements_patch_dim(geom::MappedGeometry, patch_id::Int)
    base_geom = get_base_geometry(geom, patch_id)
    return get_n_elements_patch_dim(base_geom, 1)
end

function get_n_elements_patch_dim(geom::MaskedGeometry, patch_id::Int)
    base_geom = get_base_geometry(geom)
    return get_n_elements_patch_dim(base_geom, patch_id)
end

function get_n_elements_patch_dim(geom::UnstructuredGeometry, patch_id::Int)
    patch_geom = get_geometry(geom, patch_id)
    patch_break = get_breakpoints(patch_geom)
    n_elements_dim = ntuple(d -> length(patch_break[d]) - 1, length(patch_break))
    return n_elements_dim
end 

function get_n_elements_patch_dim(
    geom::TensorProductGeometry{manifold_dim, image_dim, num_patches}, 
    patch_id::Int
    ) where {manifold_dim, image_dim, num_patches}
    constituent_geoms = get_constituent_geometries(geom)
    constituent_num_patches = map(get_num_patches, constituent_geoms)
    cart_num_patches = CartesianIndices(constituent_num_patches)
    patch_tuple = cart_num_patches[patch_id]

    dims_per_const = map(
        (g, pid) -> get_n_elements_patch_dim(g, pid),
        constituent_geoms,
        Tuple(patch_tuple)
    )
    return Tuple(vcat(dims_per_const...))
end 

function get_n_elements_patch_dim(geom::HierarchicalGeometry, level::Int, patch_id::Int)
    level_geom = get_level_geometry(geom, level)
    return get_n_elements_patch_dim(level_geom, patch_id)
end 

# function get_n_elements_patch_dim(geom::QBoxGeometry, level::Int, patch_id::Int)
#     hier_geom = get_hierarchical_geometry(geom)
#     return get_n_elements_patch_dim(hier_geom, level, patch_id)
# end 
