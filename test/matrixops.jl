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

    m = RobustModel()
    @defVar(m, matvar[1:3,1:3])
    @defUnc(m, vecunc[1:3])
    @defUnc(m, matunc[1:3,1:3])

    b = [1,2,3]
    A = eye(3,3)
    c = @addConstraint(m, A*vecunc .== b)
    for i in 1:3
        @fact conToStr(c[i]) --> "vecunc[$i] $eq $i"
    end

    c = @addConstraint(m, matunc*ones(3) .== b)
    for i in 1:3
        @fact conToStr(c[i]) --> "matunc[$i,1] + matunc[$i,2] + matunc[$i,3] $eq $i"
    end

    c = @addConstraint(m, matvar*vecunc .== b)
    for i in 1:3
        @fact conToStr(c[i]) --> "vecunc[1] matvar[$i,1] + vecunc[2] matvar[$i,2] + vecunc[3] matvar[$i,3] $eq $i"
    end

    c = @addConstraint(m, matvar*matunc .== ones(3,3))
    for i in 1:3, j in 1:3
        @fact conToStr(c[i,j]) --> "matunc[1,$j] matvar[$i,1] + matunc[2,$j] matvar[$i,2] + matunc[3,$j] matvar[$i,3] $eq 1"
    end

    c = @addConstraint(m, matvar.*matunc .== ones(3,3))
    for i in 1:3, j in 1:3
        @fact conToStr(c[i,j]) --> "matunc[$i,$j] matvar[$i,$j] $eq 1"
    end

    c = @addConstraint(m, 2.*matvar.*matunc + matvar.*matunc .== ones(3,3))
    for i in 1:3, j in 1:3
        @fact conToStr(c[i,j]) --> "(2 matunc[$i,$j]) matvar[$i,$j] + matunc[$i,$j] matvar[$i,$j] $eq 1"
    end

end