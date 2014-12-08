using JuMPeR
using Base.Test

# Use Gurobi it at all possible
solver = JuMPeR.JuMP.UnsetSolver()
if Pkg.installed("Gurobi") != nothing
    using Gurobi
    solver = GurobiSolver(OutputFlag=0)
    println("Selected Gurobi as solver")
end

tests = [   "operators.jl",
            #"macro.jl",
            #"polyhedral.jl",
            #"bertsim.jl",
            #"oracle.jl",
            #"ellipse.jl",
            #"scenario.jl",
            #"resolve.jl",
            #"graph.jl"]
        ]

println("Running tests...")
for curtest in tests
    println("Test: $curtest")
    include(curtest)
end