# Adaptive Loop

This folder demonstrates adaptive mesh refinement.
It solves Poisson's equation on a 2D Cartesian domain using hierarchical B-spline spaces,
refining the mesh adaptively based on the error at each step.

## Running the adaptive loop

From this directory, run the script using Julia with the local project environment
activated:

```bash
julia --project adaptive-loop.jl
```

Or, from the Julia REPL:

```julia
# Activate the environment first
import Pkg; Pkg.activate(".")

# Then run the script
include("adaptive-loop.jl")
```

## Output

The script writes the computed and analytical solutions to VTK files in the current
directory:

- `Adaptive-Poisson-uₕ.vtu` — computed solution
- `Adaptive-Poisson-u.vtu` — analytical solution

These can be visualised with [ParaView](https://www.paraview.org/), for example.
