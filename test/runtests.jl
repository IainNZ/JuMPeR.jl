using JuMPeR
using Base.Test
solver = nothing
if Pkg.installed("Gurobi") != nothing
    using Gurobi
    solver = GurobiSolver(OutputFlag=0)
    println("Selected Gurobi as solver")
end

tests = [   "operators.jl",
            "polyhedral.jl",
            "bertsim.jl",
            "affine.jl",
            "oracle.jl",
            "ellipse.jl",
            "scenario.jl",
            "graph.jl"]

println("Running tests...")
for curtest in tests
    println("Test: $curtest")
    include(curtest)
end