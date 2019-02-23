#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2016: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# test/adp_newsvendor.jl
# Tests adaptive code for a two-stage newsvendor-style problem.
# We have demand at N locations (Dᵢ), and much satisfy demand by buying
# inventory now (x). Amount sold is then no more than demand, and total
# amount sold must be less than amount bought initially.
#
# max z
# st  z ≤ ∑S(D) - 0.5x    ∀D
#     x ≥ ∑S(D)
#     0 ≤ Sᵢ(D) ≤ Dᵢ
#     D ∈ { ‖(D - μ)/σ‖₁ ≤ ⌊√N⌋, ‖(D - μ)/σ‖∞ ≤ 1 }
#-----------------------------------------------------------------------

using JuMP, JuMPeR
using Test

const TOL = 1e-4

if !(:lp_solvers in names(Main))
    printstyled("Loading solvers...\n", color = :magenta)
    include(joinpath(dirname(pathof(JuMP)),"..","test","solvers.jl"))
end
lp_solvers  = filter(s->(!occursin("SCSSolver", string(typeof(s)))), lp_solvers)
solver_name(solver) = split(string(typeof(solver)),".")[2]

@testset "Adaptive Newsvendor Model" begin
printstyled("Adaptive Newsvendor Model...\n", color = :yellow)
@testset "with $(solver_name(solver)), cuts=$cuts" for
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
        @uncertain(m, D[1:N])
        @constraint(m, norm((D - Dμ)./Dσ,   1) <= isqrt(N))
        @constraint(m, norm((D - Dμ)./Dσ, Inf) <= 1)
        return D
    end

    @testset "Static, manual" begin
        m = RobustModel(solver=solver)
        D = add_set(m)
        @variable(m, x >= 0)
        @variable(m, S[1:N] >= 0)
        @variable(m, z <= 1000);  @objective(m, Max, z)
        @constraint(m, z <= sum(S) - 0.5*x)
        for i in 1:N; @constraint(m, S[i] <= D[i]); end
        @constraint(m, sum(S) <= x)
        solve(m, prefer_cuts=cuts)
        @test isapprox(getvalue(x), static_x, atol=TOL)
        @test isapprox(getvalue(z), static_z, atol=TOL)
    end  # "Static, manual"

    @testset "Static, auto" begin
        m = RobustModel(solver=solver)
        D = add_set(m)
        @variable(m, x >= 0)
        @adaptive(m, S[1:N] >= 0, policy=Static, depends_on=D)
        @variable(m, z <= 1000);  @objective(m, Max, z)
        @constraint(m, z <= sum(S) - 0.5*x)
        for i in 1:N; @constraint(m, S[i] <= D[i]); end
        @constraint(m, sum(S) <= x)
        solve(m, prefer_cuts=cuts)
        @test isapprox(getvalue(x), static_x, atol=TOL)
        @test isapprox(getvalue(z), static_z, atol=TOL)
    end  # "Static, auto"

    @testset "Affine, manual" begin
        m = RobustModel(solver=solver)
        D = add_set(m)
        @variable(m, x >= 0)
        @variable(m, S_aff[1:N,0:N])
        S = Any[S_aff[i,0] for i in 1:N]
        for i in 1:N
            for j in 1:N
                S[i] += S_aff[i,j] * D[j]
            end
            @constraint(m, S[i] >= 0)
        end
        @variable(m, z <= 1000);  @objective(m, Max, z)
        @constraint(m, z <= sum(S) - 0.5*x)
        for i in 1:N; @constraint(m, S[i] <= D[i]); end
        @constraint(m, sum(S) <= x)
        solve(m, prefer_cuts=cuts)
        @test isapprox(getvalue(x), aff_x, atol=TOL)
        @test isapprox(getvalue(z), aff_z, atol=TOL)
    end

    @testset "Affine, auto" begin
        m = RobustModel(solver=solver)
        D = add_set(m)
        @variable(m, x >= 0)
        @adaptive(m, S[1:N] >= 0, policy=Affine, depends_on=D)
        @variable(m, z <= 1000);  @objective(m, Max, z)
        @constraint(m, z <= sum(S) - 0.5*x)
        for i in 1:N; @constraint(m, S[i] <= D[i]); end
        @constraint(m, sum(S) <= x)
        solve(m, prefer_cuts=cuts)
        @test isapprox(getvalue(x), aff_x, atol=TOL)
        @test isapprox(getvalue(z), aff_z, atol=TOL)
    end

end  # "with ..."
end  # "Newsvendor"
