using JuMP, JuMPeR
using FactCheck
using BaseTestNext

# Create list of solvers using JuMP's code
print_with_color(:yelow, "Loading solvers...\n")
include(joinpath(Pkg.dir("JuMP"),"test","solvers.jl"))

@testset "JuMPeR" begin
    include("operators.jl")
end

tests=[ #"operators.jl",
        #"matrixops.jl",
        "print.jl",
        "macro.jl",
        "oracle.jl",
        "oracle_general.jl",
        "oracle_general_L1.jl",
        "oracle_general_Linf.jl",
        "oracle_bertsim.jl",
        "oracle_general_graph.jl",
        "scenario.jl"]

println("Running tests...")
for curtest in tests
    include(curtest)
end

FactCheck.exitstatus()
