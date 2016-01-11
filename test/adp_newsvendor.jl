#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2015: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# test/adp_newsvendor.jl
# Tests adaptive code on the two-stage newsvendor problem
# max z
# st  z ≤ ∑S(D) - 0.5x    ∀D
#     x ≥ ∑S(D)
#     0 ≤ Sᵢ(D) ≤ Dᵢ
#     D ∈ { ‖(D - μ)/σ‖₁ ≤ ⌊√N⌋, ‖(D - μ)/σ‖∞ ≤ 1 }
#-----------------------------------------------------------------------

using JuMP, JuMPeR
using BaseTestNext

const TOL = 1e-4

if !(:lp_solvers in names(Main))
    print_with_color(:yellow, "Loading solvers...\n")
    include(joinpath(Pkg.dir("JuMP"),"test","solvers.jl"))
end
lp_solvers  = filter(s->(!contains(string(typeof(s)),"SCSSolver")), lp_solvers)


@testset "Newsvendor" begin
print_with_color(:yellow, "Newsvendor...\n")

@testset "With lp_solver $(typeof(solver)), cuts=$cuts" for
                        solver in lp_solvers, cuts in [true,false]

    # Generate data
    N = 5
    Dμ = [28,32,38,40,46]
    Dσ = 0.5 * Dμ

    # Calculate solutions manually
    static_x = sum(Dμ-Dσ)
    static_z = 0.5*static_x
    aff_x = sum(Dμ[1:N-isqrt(N)]) + 0.5*sum(Dμ[N-isqrt(N)+1:N])
    aff_z = 0.5*aff_x

    function add_set(m)
        @defUnc(m, D[1:N])
        @addConstraint(m, norm((D - Dμ)./Dσ,   1) <= isqrt(N))
        @addConstraint(m, norm((D - Dμ)./Dσ, Inf) <= 1)
        return D
    end

    @testset "Static policy, manual" begin
        m = RobustModel()
        D = add_set(m)
        @defVar(m, x >= 0)
        @defVar(m, S[1:N] >= 0)
        @defVar(m, z <= 1000);  @setObjective(m, Max, z)
        @addConstraint(m, z <= sum(S) - 0.5*x)
        for i in 1:N; @addConstraint(m, S[i] <= D[i]); end
        @addConstraint(m, sum(S) <= x)
        solve(m, prefer_cuts=cuts)
        @test isapprox(getValue(x), static_x, atol=TOL)
        @test isapprox(getValue(z), static_z, atol=TOL)
    end  # "Static, manual"

    @testset "Static, auto" begin
        m = RobustModel()
        D = add_set(m)
        @defVar(m, x >= 0)
        @defAdaptVar(m, S[1:N] >= 0, policy=Static, depends_on=D)
        @defVar(m, z <= 1000);  @setObjective(m, Max, z)
        @addConstraint(m, z <= sum(S) - 0.5*x)
        for i in 1:N; @addConstraint(m, S[i] <= D[i]); end
        @addConstraint(m, sum(S) <= x)
        solve(m, prefer_cuts=cuts)
        @test isapprox(getValue(x), static_x, atol=TOL)
        @test isapprox(getValue(z), static_z, atol=TOL)
    end  # "Static, auto"

    @testset "Affine, manual" begin
        m = RobustModel()
        D = add_set(m)
        @defVar(m, x >= 0)
        @defVar(m, S_aff[1:N,0:N])
        S = Any[S_aff[i,0] for i in 1:N]
        for i in 1:N
            for j in 1:N
                S[i] += S_aff[i,j] * D[j]
            end
            @addConstraint(m, S[i] >= 0)
        end
        @defVar(m, z <= 1000);  @setObjective(m, Max, z)
        @addConstraint(m, z <= sum(S) - 0.5*x)
        for i in 1:N; @addConstraint(m, S[i] <= D[i]); end
        @addConstraint(m, sum(S) <= x)
        solve(m, prefer_cuts=cuts)
        @test isapprox(getValue(x), aff_x, atol=TOL)
        @test isapprox(getValue(z), aff_z, atol=TOL)
    end

    @testset "Affine, auto" begin
        m = RobustModel()
        D = add_set(m)
        @defVar(m, x >= 0)
        @defAdaptVar(m, S[1:N] >= 0, policy=Affine, depends_on=D)
        @defVar(m, z <= 1000);  @setObjective(m, Max, z)
        @addConstraint(m, z <= sum(S) - 0.5*x)
        for i in 1:N; @addConstraint(m, S[i] <= D[i]); end
        @addConstraint(m, sum(S) <= x)
        solve(m, prefer_cuts=cuts)
        @test isapprox(getValue(x), aff_x, atol=TOL)
        @test isapprox(getValue(z), aff_z, atol=TOL)
    end

end  # "lp_solver"
end  # "Newsvendor"
