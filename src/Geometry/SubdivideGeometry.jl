"""
	subdivide_geometry(parent_geo::AbstractGeometry, num_subdivisions)

Returns a refined version of `parent_geo` where each element is subdivided according to
`num_subdivisions`.

# Returns
- `child_geometry`: The geometry corresponding to the subdivision of `parent_geo`.
"""
function subdivide_geometry(parent_geo::AbstractGeometry, num_subdivisions)
    throw(MethodError(subdivide_geometry, (parent_geo, num_subdivisions)))
end

function subdivide_geometry(
    parent_geo::UnstructuredGeometry{manifold_dim}, num_subdivisions::Int
) where {manifold_dim}
    return subdivide_geometry(parent_geo, ntuple(dim -> num_subdivisions, manifold_dim))
end

function subdivide_geometry(
    parent_geo::UnstructuredGeometry{manifold_dim, image_dim, num_patches},
    num_subdivisions::NTuple{manifold_dim},
) where {manifold_dim, image_dim, num_patches}
    return subdivide_geometry(parent_geo, ntuple(patch -> num_subdivisions, num_patches))
end

function subdivide_geometry(
    parent_geo::UnstructuredGeometry{manifold_dim, image_dim, num_patches},
    num_subdivisions::NTuple{num_patches, NTuple{manifold_dim}},
) where {manifold_dim, image_dim, num_patches}
    patch_parent_geo = get_geometry_per_patch(parent_geo)
    patch_child_geo = ntuple(
        patch -> subdivide_geometry(patch_parent_geo[patch], num_subdivisions[patch]),
        num_patches,
    )
    child_geo = UnstructuredGeometry(patch_child_geo)

    return child_geo
end

function subdivide_geometry(hier_geo::HierarchicalGeometry, num_subdivisons)
    parent_geo = get_level_geometry(hier_geo, get_num_levels(hier_geo))
    child_geo = subdivide_geometry(parent_geo, num_subdivisons)
    new_geometries = (get_geometries(hier_geo)..., child_geo)
    new_active = (get_active_elements(hier_geo)..., Int[])
    new_hier_geo = HierarchicalGeometry(new_geometries, Hierarchy.ActiveInfo(new_active))
    
    return new_hier_geo
end

function subdivide_geometry!(qbox_geo::QBoxGeometry, num_subdivisons)
    hier_geo = get_hierarchical_geometry(qbox_geo)
    new_hier_geo = subdivide_geometry(hier_geo, num_subdivisons)
    qbox_geo.hier_geom = new_hier_geo
    #get_hierarchical_geometry(qbox_geo) = new_hier_geo
end

function subdivide_geometry(parent_geo::MappedGeometry, num_subdivisons)
    parent_base_geometry = get_base_geometry(parent_geo)
    child_base_geometry = subdivide_geometry(parent_base_geometry, num_subdivisons)
    mapping = get_mapping(parent_geo)
    child_geometry = MappedGeometry(child_base_geometry, mapping)

    return child_geometry
end

function subdivide_geometry(parent_geo::MaskedGeometry, num_subdivisons)
    parent_base_geometry = get_base_geometry(parent_geo)
    child_base_geometry = subdivide_geometry(parent_base_geometry, num_subdivisons)
    mask = get_evaluation_mask(parent_geo)
    child_geometry = MaskedGeometry(child_base_geometry, mask)

    return child_geometry
end

function subdivide_geometry(
    parent_geo::TensorProductGeometry{
        manifold_dim, image_dim, num_patches, num_geometries
    },
    num_subdivisons::NTuple{manifold_dim},
) where {manifold_dim, image_dim, num_patches, num_geometries}
    const_parent_geo = get_constituent_geometries(parent_geo)
    const_manifold_indices = get_constituent_manifold_indices(parent_geo)
    const_num_subdivions = ntuple(
        geo -> num_subdivisons[const_manifold_indices[geo]], num_geometries
    )
    const_child_geo = ntuple(
        geo -> subdivide_geometry(const_parent_geo[geo], const_num_subdivions[geo]),
        num_geometries,
    )
    child_geometry = TensorProductGeometry(const_child_geo)

    return child_geometry
end

function subdivide_geometry(parent_geo::CartesianGeometry{1}, num_subdivisons::Int)
    return subdivide_geometry(parent_geo, (num_subdivisons,))
end

function subdivide_geometry(
    parent_geo::CartesianGeometry{manifold_dim, image_dim, 1},
    num_subdivisons::NTuple{manifold_dim},
) where {manifold_dim, image_dim}
    parent_breakpoints = get_breakpoints(parent_geo)
    const_child_breakpoints = ntuple(
        dim -> subdivide_breakpoints(parent_breakpoints[dim], num_subdivisons[dim]),
        manifold_dim,
    )
    child_geometry = CartesianGeometry(const_child_breakpoints)

    return child_geometry
end

function subdivide_geometry(
    parent_geo::CartesianGeometry{manifold_dim, image_dim, num_patches},
    num_subdivisions::NTuple{manifold_dim}
) where {manifold_dim, image_dim, num_patches}
    return subdivide_geometry(parent_geo, ntuple(patch -> num_subdivisions, num_patches))
end

function subdivide_geometry(
    parent_geo::CartesianGeometry{manifold_dim, image_dim, num_patches},
    num_subdivisions::NTuple{num_patches, NTuple{manifold_dim}}
) where {manifold_dim, image_dim, num_patches}

    child_breaks = ntuple(patch -> begin
        patch_breaks = get_breakpoints(parent_geo, patch)
        patch_subdivs = num_subdivisions[patch]

        ntuple(dim -> subdivide_breakpoints(
            patch_breaks[dim],
            patch_subdivs[dim]
        ), manifold_dim)

    end, num_patches)

    return CartesianGeometry(child_breaks)
end

"""
    subdivide_breakpoints(parent_breakpoints::Vector{Float64}, num_subdivisions::Int)

Subdivides `parent_breakpoints` by uniformly subdiving each element `num_subdivisions`
times.

# Arguments
- `parent_breakpoints::AbstractVector`: Parent set of breakpoints.
- `num_subdivisions`: Number of times each element is subdivided.
# Returns
- `child_breakpoints::Vector{Float64}`: Child set of breakpoints.
"""
function subdivide_breakpoints(parent_breakpoints::AbstractVector, num_subdivisions)
    num_parent_breakpoints = length(parent_breakpoints)
    num_child_breakpoints =
        num_parent_breakpoints + (num_parent_breakpoints - 1) * (num_subdivisions - 1)
    child_breakpoints = Vector{Float64}(undef, num_child_breakpoints)
    child_breakpoints[end] = parent_breakpoints[end]
    step_size = 1 / num_subdivisions
    for i in 1:(num_parent_breakpoints - 1), j in 0:(num_subdivisions - 1)
        index = (i - 1) * num_subdivisions + j + 1
        child_breakpoints[index] =
            parent_breakpoints[i] +
            j * step_size * (parent_breakpoints[i + 1] - parent_breakpoints[i])
    end

    return child_breakpoints
end

function subdivide_breakpoints(parent_breakpoints::LinRange, num_subdivisions)
    start_breakpoint = first(parent_breakpoints)
    end_breakpoint = last(parent_breakpoints)
    parent_num_elements = length(parent_breakpoints) - 1
    child_breakpoints = LinRange(
        start_breakpoint, end_breakpoint, parent_num_elements * num_subdivisions + 1
    )

    return child_breakpoints
end
