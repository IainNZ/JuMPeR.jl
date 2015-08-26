#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2015: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# test/matrixops.jl
# Test doing matrix-variable-uncertain operations
#-----------------------------------------------------------------------

using JuMP, JuMPeR
using FactCheck

# To ensure the tests work on both Windows and Linux/OSX, we need
# to use the correct comparison operators in the strings
const leq = JuMP.repl[:leq]
const geq = JuMP.repl[:geq]
const  eq = JuMP.repl[:eq]

facts("[matrixops] Matrix operation tests") do

context("Numbers with Uncertains") do
    m = RobustModel()
    @defUnc(m, vecunc[1:3])
    @defUnc(m, matunc[1:3,1:3])
    b = rand(3)
    A = rand(3,3)
    @addConstraint(m, 0 .== A*vecunc + b)
end



end