module QBoxGeometryRefinementTests

using Mantis
using Test

############################################################################################
#                                          Setup                                           #
############################################################################################

#2D QBox
geom_2D = Geometry.CartesianGeometry((
    0.0:2.0:5.0,
    0.0:2.0:5.0
))
n1_2D = Geometry.get_num_elements(geom_2D)   # = 16
qbox_size_2D = (2,2)
num_subdivisions_2D = (2,1)
qbg_2D = Geometry.QBoxGeometry_refine(geom_2D, qbox_size_2D, num_subdivisions_2D)
#Plot.plot(qbg_2D; vtk_filename="Starting Geometry Refinement Tests")

@testset "QBox Refinement tests 2D" begin
    error = [0.05, 0.05, 0.4, 0.3, 0.05, 0.05, 0.3, 0.4, 0.05, 0.8, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05]
    remove = Geometry.refine_qboxgeom_max!(qbg_2D, error, 0.1 )
    #@test remove == Set([(3,1,1)])
    @test remove == [[9,10,13,14]]
    #Plot.plot(qbg_2D; vtk_filename="Refinement after max")
    error2 = [0.05, 0.05, 0.4, 0.3, 0.05, 0.05, 0.3, 0.4,0.05, 0.05, 0.05, 0.05, 0.35, 0.33,0.34, 0.34, 0.05, 0.05, 0.05, 0.05]
    remove2 = Geometry.refine_qboxgeom_avg!(qbg_2D, error2, 0.1 )
    #@test remove2 == [(5,2,1),(2,1,1)]
    @test remove2 == [[3,4,7,8],[17,18,25,26]]
    #Plot.plot(qbg_2D; vtk_filename="Refinement after avg")
end

end