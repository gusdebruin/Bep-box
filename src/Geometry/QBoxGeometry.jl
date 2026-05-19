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
                    "Number of 'active_elements' must be divisible by the qbox_size." 
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

function get_hierarchical_geometry(geometry::QBoxGeometry)
    return geometry.hier_geom
end

function get_qbox_size(geometry::QBoxGeometry)
    return geometry.qbox_size
end
