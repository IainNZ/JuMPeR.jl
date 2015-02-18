using JuMP, JuMPeR
using Base.Test, FactCheck
FactCheck.setstyle(:compact)

# Create list of solvers using JuMP's code
println("Loading solvers...")
include(joinpath(Pkg.dir("JuMP"),"test","solvers.jl"))

tests=[ "operators.jl",
        "print.jl",
        "macro.jl",
        "oracle.jl",
        "oracle_general.jl",
        "oracle_bertsim.jl",
        "oracle_general_graph.jl",
        "scenario.jl"]

println("Running tests...")
for curtest in tests
    include(curtest)
end

FactCheck.exitstatus()