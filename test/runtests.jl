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
# Launch all JuMPeR tests - run with Pkg.test("JuMPeR")
#-----------------------------------------------------------------------

using JuMP, JuMPeR
using Test, LinearAlgebra

# Use JuMP's testing code to load available solvers
# and provide vectors of solvers by capability
printstyled("Loading solvers...\n", color = :magenta)
include(joinpath(dirname(pathof(JuMP)),"..","test","solvers.jl"))

@testset "JuMPeR" begin
    include("operators.jl")
    include("print.jl")
    include("macro.jl")
    include("uncsets.jl")
    include("uncsets_basic.jl")
    include("uncsets_basic_L1.jl")
    include("uncsets_basic_L2.jl")
    include("uncsets_basic_Linf.jl")
    include("uncsets_budget.jl")
    include("adp_newsvendor.jl")
    include("adp_inventory.jl")
end
