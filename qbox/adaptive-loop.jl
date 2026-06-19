using Mantis
import CairoMakie as CM
############################################################################################
#                                   Problem Description                                    #
############################################################################################

#=

We want to solve Poisson's equation using adaptive refinement. The weak formulation of the
problem can be written as: Find u ∈ ℍ such that for f ∈ L²

    ∫ ∇v⋅∇u dΩ = ∫ v⋅f dΩ, ∀ v ∈ ℍ. 

We start by solving the problem on an initial, coarse level, and then perform refinement
based on the error of the computed solution. The overall idea looks like this:

    Solve -> Estimate Errors -> Mark for Refinement -> Refine space -> Solve -> ...

In this case, we will know that the analytical solution is, so the `Estimate Errors` step is
actually a `Compute Errors` step.


When working with q-boxes, the `Mark for Refinement` step is the interesting one: that is
where q-boxes come in to define what, and how, elements can be refined.

=#

############################################################################################
#                                      Problem setup                                       #
############################################################################################

#= 
We start by defining the geometry for the first level. In this case, we want a 2D Cartesian
geometry:
=#
const starting_point = (0.0, 0.0)
const box_size = (1.0, 1.0)
const num_elements = (12, 12)
#=
Then, we define the parameters for the B-spline spaces of each level:
=#
const p = (2, 2) # polynomial degree
const k = (1, 1) # regularity
#=
Because we will want to adaptively refine in a loop, we create the B-spline space and
immediately turn it into a trivial hierarchical space with no refinement. For that we just
need to specify the number of subdivisions that will be used:
=#
num_subdivisions = (2, 2)
B = FunctionSpaces.create_bspline_space(starting_point, box_size, num_elements, p, k)
H = FunctionSpaces.HierarchicalFiniteElementSpace(B, num_subdivisions)

Mantis.Plot.plot(FunctionSpaces.get_geometry(H); vtk_filename="Starting Geometry Poisson no QBox")

#=
To set up our adaptive loop we also need to define a few other things:

1. The number of steps in the adaptive loop.
2. The Dörfler parameter to specify how many elements get refined.
=#
N = 10 # Number of steps.
θ = 0.20 # Dörfler parameter.

############################################################################################
#                                          Solve                                           #
############################################################################################

#=
Before starting our loop, we of course need to define what the forcing term f is. The
expression we will use looks kind of complicated, so that we get an interesting refinement
pattern. Because the details of how the function is defined are not very important, we will
put it in the file `forcing.jl` and just load it here.

We just need to know that it gives as a function `forcing` that we can call with a
`geometry` to obtain our forcing term: `forcing(geometry)`.

For similar reasons, we handle the analytical solution the same way, by loading it from the
file `solution.jl`.

Both the forcing, and by consequence the solution, will depend on a parameter ϵ that
determines how “sharp” the functions are — by sharp here, we mean how local the features
are. In this case, the lower the value of ϵ the “sharper” the function. It's expected that
local refinement will concentrate around the sharp feature.
=#
const ϵ = 0.01
include("forcing.jl")
include("solution.jl")
#= 
For example, we could already take a look at the analytical solution by running:
```julia
U = solution(FunctionSpaces.get_geometry(B))
Plot.export_form_fields_to_vtk((U,), "solution")
```
=#

#=
The last step before writing the adaptive loop is just defining quadrature. This is what
tells `Mantis` what points should be used to compute the integrals in the weak formulation.
=#
qr = first(Quadrature.get_canonical_quadrature_rules(Quadrature.gauss_legendre, p .+ 1))

#=
We are now finally ready to write the adaptive loop! The function does a lot of things, so
we will include various comments inside the function.
=#
function adaptive_loop(H, N, θ)
    # To keep track of what our code is doing, we will print things like the line below:
    println("Solving the problem on the initial step...")
    dofs_history = Int[]
    error_history = Float64[]

    #=
    Before we defined `qr`, which is the quadrature rule for a single element. Now we just
    need to say what we will used that same rule for every element in our geometry:
    =#
    geometry = FunctionSpaces.get_geometry(H) # Get the current geometry.
    Q = Quadrature.StandardQuadrature(qr, Geometry.get_num_elements(geometry))

    # Next, we define the forcing and analytical solution on the current geometry
    f = forcing(geometry)
    u = solution(geometry)

    #=
    Then, we solve the problem at the current step. Because `Mantis` only works with
    differential forms, we need to wrap or function space `H` as a `FormSpace`. This is a
    minor detail, we can still think of this `FormSpace` as our function space.
    =#
    U = Forms.FormSpace(0, H, "uₕ")
    #=
    This will solve Poisson's equation, as we described above; the term “zero form Hodge
    Laplacian” is just another name for the same the equation.
    =#
    uₕ = Assemblers.solve_zero_form_hodge_laplacian(U, f, Q)

    #=
    Finally, because we know the analytical solution, computing the error per element is
    easy.
    =#
    err = Analysis.compute_error_per_element(uₕ, u, Q)
    # We can also take a look at what is the total error we got:
    total_err = sqrt(sum(e -> e^2, err))
    push!(dofs_history, FunctionSpaces.get_num_basis(H))
    push!(error_history, total_err)
    @show total_err

    # Now we just repeat what we did above, while refining the space at each step.
    for step in 1:N
        println("Solving the problem on step $(step)...")
        #=
        This is a specific pattern of refinement. It will refine the support of all basis
        functions whose supports intersect the `dorfler_marking`.
        Q-boxes will likely do something different.
        =#
        dorfler_marking = FunctionSpaces.get_dorfler_marking(err, θ)
        marked_elements_per_level = FunctionSpaces.get_padding_per_level(
            Forms.get_fe_space(U), dorfler_marking
        )

        # And we use the selected elements to refine the hierarchical space.
        H = FunctionSpaces.refine_space(
            Forms.get_fe_space(U), marked_elements_per_level
        )

        # Then, we repeat what we did at the initial step
        geometry = FunctionSpaces.get_geometry(H) # Get the current geometry.
        Q = Quadrature.StandardQuadrature(qr, Geometry.get_num_elements(geometry))
        f = forcing(geometry)
        u = solution(geometry)

        U = Forms.FormSpace(0, H, "uₕ")
        uₕ = Assemblers.solve_zero_form_hodge_laplacian(U, f, Q)

        err = Analysis.compute_error_per_element(uₕ, u, Q)
        total_err = sqrt(sum(e -> e^2, err))
        push!(dofs_history, FunctionSpaces.get_num_basis(H))
        push!(error_history, total_err)
        @show total_err
    end
    @show dofs_history
    @show error_history

    return uₕ, u, dofs_history, error_history
end

############################################################################################
#                                     Getting Results                                      #
############################################################################################


#=
To get the results from our adaptive loop we just need to call it. It will give us both the
computed and analytical solutions.
=#
uₕ, u, dofs_history, error_history = adaptive_loop(H, N, θ)

# To finish it off we export the results so we can take a look at them.
Plot.export_form_fields_to_vtk((uₕ, u), "Adaptive-Poisson no QBox")


# Plot 1: DOFs vs fout (log-log)
fig1 = CM.Figure()
ax1 = CM.Axis(fig1[1,1],
    xlabel = "Total DOFs",
    ylabel = "L² error",
    title  = "Convergence adaptive refinement",
    yscale = log10
)
CM.lines!(ax1, dofs_history, error_history)
CM.scatter!(ax1, dofs_history, error_history)
CM.save("convergentie no QBox.png", fig1)