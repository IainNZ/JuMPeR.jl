#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2016: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# test/print.jl
# Test pretty-printing-related functionality not otherwise covered
#-----------------------------------------------------------------------

using JuMP, JuMPeR
using BaseTestNext
import JuMP: repl

# Helper function to test IO methods work correctly
function io_test(obj, exp_str; repl=:both)
    @test sprint(print, obj) == exp_str
end

@testset "Printing" begin
print_with_color(:yellow, "Printing...\n")

@testset "RobustModel" begin
    le, ge, fa = repl[:leq], repl[:geq], repl[:for_all]
    inset, dots = repl[:in], repl[:dots]
    infty, union = repl[:infty], repl[:union]

    mod_1 = RobustModel()
    @variable(mod_1, vars[1:10])

    @uncertain(mod_1, a>=1)
    @uncertain(mod_1, b<=1)
    @uncertain(mod_1, -1<=c<=1)
    @uncertain(mod_1, a1>=1,Int)
    @uncertain(mod_1, b1<=1,Int)
    @uncertain(mod_1, -1<=c1<=1,Int)
    @uncertain(mod_1, x, Bin)
    @uncertain(mod_1, y)
    @uncertain(mod_1, z, Int)
    @uncertain(mod_1,      bnd_free[2:5])
    @uncertain(mod_1,      bnd_lowb[2:5] >= 2)
    @uncertain(mod_1,      bnd_high[2:5] <= 5)
    @uncertain(mod_1, 2 <= bnd_both[2:5] <= 5)
    @uncertain(mod_1, mat2d[1:3,1:3])
    @uncertain(mod_1, mat3d[1:3,1:3,1:3])

    @objective(mod_1, Max, 2*vars[1])
    # Deterministic
    @constraint(mod_1, vars[10] <= 10)
    # Mixed
    @constraint(mod_1, a*vars[5] <= 5)
    # Uncertain
    # @constraint(mod_1, a + b <= 2)
    # a + b $le 2
    # Ellipse
    # @constraint(mod_1, norm([a,b]) <= 1)
    # ‖a,b‖₂ $le 1

    io_test(mod_1, """
Max 2 vars[1]
Subject to
 vars[10] ≤ 10
 vars[i] ∀ i ∈ {1,2,…,9,10}
Uncertain constraints:
a vars[5] ≤ 5
Uncertain parameters:
bnd_free[i] free ∀ i ∈ {2,3,4,5}
bnd_lowb[i] ≥ 2 ∀ i ∈ {2,3,4,5}
bnd_high[i] ≤ 5 ∀ i ∈ {2,3,4,5}
2 ≤ bnd_both[i] ≤ 5 ∀ i ∈ {2,3,4,5}
mat2d[i,j] free ∀ i ∈ {1,2,3}, j ∈ {1,2,3}
mat3d[i,j,k] free ∀ i ∈ {1,2,3}, j ∈ {1,2,3}, k ∈ {1,2,3}
a ≥ 1
b ≤ 1
-1 ≤ c ≤ 1
a1 ≥ 1, integer
b1 ≤ 1, integer
-1 ≤ c1 ≤ 1, integer
x ∈ {0,1}
y free
z free, integer
""", repl=:print)
end  # "RobustModel"

@testset "JuMPContainer{Uncertain}" begin
    le, ge, fa = repl[:leq], repl[:geq], repl[:for_all]
    inset, dots = repl[:in], repl[:dots]
    infty, union = repl[:infty], repl[:union]

    m = RobustModel()
    @uncertain(m,      bnd_free[2:5])
    @uncertain(m,      bnd_lowb[2:5] >= 2)
    @uncertain(m,      bnd_high[2:5] <= 5)
    @uncertain(m, 2 <= bnd_both[2:5] <= 5)
    @uncertain(m,      bnd_difflo[i=2:5] >= i)
    @uncertain(m,      bnd_diffup[i=2:5] <= i)
    @uncertain(m, i <= bnd_diffbo[i=2:5] <= 2i)
    @uncertain(m, i <= bnd_difflo_with_up[i=2:5] <= 5)
    @uncertain(m, 2 <= bnd_diffup_with_lo[i=2:5] <= i)

    io_test(bnd_free, "bnd_free[i] free $fa i $inset {2,3,4,5}")
    io_test(bnd_lowb, "bnd_lowb[i] $ge 2 $fa i $inset {2,3,4,5}")
    io_test(bnd_high, "bnd_high[i] $le 5 $fa i $inset {2,3,4,5}")
    io_test(bnd_both, "2 $le bnd_both[i] $le 5 $fa i $inset {2,3,4,5}")
    io_test(bnd_difflo, "bnd_difflo[i] $ge $dots $fa i $inset {2,3,4,5}")
    io_test(bnd_diffup, "bnd_diffup[i] $le $dots $fa i $inset {2,3,4,5}")
    io_test(bnd_diffbo, "$dots $le bnd_diffbo[i] $le $dots $fa i $inset {2,3,4,5}")
    io_test(bnd_difflo_with_up, "$dots $le bnd_difflo_with_up[i] $le 5 $fa i $inset {2,3,4,5}")
    io_test(bnd_diffup_with_lo, "2 $le bnd_diffup_with_lo[i] $le $dots $fa i $inset {2,3,4,5}")
end  # "JuMPContainer{Uncertain}"

end  # "Printing"
