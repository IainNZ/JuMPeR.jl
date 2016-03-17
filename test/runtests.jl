#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2016: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# test/runtests.jl
# Launch all JuMPeR tests. Run by Pkg.test("JuMPeR")
#-----------------------------------------------------------------------

using JuMP, JuMPeR
using FactCheck
using BaseTestNext

# Use JuMP's testing code to load available solvers
# and provide vectors of solvers by capability
print_with_color(:magenta, "Loading solvers...\n")
include(joinpath(Pkg.dir("JuMP"),"test","solvers.jl"))

@testset "JuMPeR" begin
    include("operators.jl")
    include("print.jl")
    include("macro.jl")
    include("oracle.jl")
    include("oracle_general.jl")
end

tests=[ "oracle_general_L1.jl",
        "oracle_general_L2.jl",
        "oracle_general_Linf.jl",
        "oracle_bertsim.jl",
        "oracle_general_graph.jl",
        "scenario.jl"]

println("Running tests...")
for curtest in tests
    include(curtest)
end

FactCheck.exitstatus()
