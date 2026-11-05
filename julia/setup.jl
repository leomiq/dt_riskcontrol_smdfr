# Instantiate the Julia environment for the monolithic MPC formulation (F1-F3).
# Usage, from the repository root:  julia --project=julia julia/setup.jl
using Pkg
Pkg.add(["JuMP", "Ipopt", "ForwardDiff", "MAT"])
Pkg.precompile()
Pkg.status()
