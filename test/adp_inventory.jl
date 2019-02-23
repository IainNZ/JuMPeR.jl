#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2016: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# test/adp_inventory.jl
# Tests adaptive code on the multi-stage inventory problem from
#   A. Ben-Tal, A. Goryashko, E. Guslitzer, A. Nemirovski
#   "Adjustable Robust Solutions of Uncertain Linear Programs"
#-----------------------------------------------------------------------

using JuMP, JuMPeR
using Test

TOL = 5e-3

if !(:lp_solvers in names(Main))
    printstyled("Loading solvers...\n", color = :magenta)
    include(joinpath(dirname(pathof(JuMP)),"..","test","solvers.jl"))
end
lp_solvers  = filter(s->(!occursin("SCSSolver", string(typeof(s)))), lp_solvers)
solver_name(solver) = split(string(typeof(solver)),".")[2]

@testset "Adaptive Inventory" begin
printstyled("Adaptive Inventory Model...\n", color = :yellow)
@testset "with $(solver_name(solver)), cuts=$cuts" for
            solver in lp_solvers, cuts in [true,false]

    #-----------------------------------------------------------------------
    # Data
    I = 3           # Number of factories
    T = 24          # Number of time periods
    d_nom = 1000*[1 + 0.5*sin(π*(t-1)/12)  for t in 1:T]  # Nominal demand
    θ = 0.20        # Uncertainty level
    α = [1.0, 1.5, 2.0]  # Production costs
    c = [α[i] * (1 + 0.5*sin(π*(t-1)/12))  for i in 1:I, t in 1:T]
    P = 567         # Maximimum production per period
    Q = 13600       # Maximumum production over all
    Vmin = 500      # Minimum inventory at warehouse
    Vmax = 2000     # Maximum inventory at warehouse
    v1 = Vmin       # Initial inventory (not provided in paper)

    function add_constraints(rm, p, d)
        # Constraint: cannot exceed production limits
        for i in 1:I, t in 1:T
            @constraint(rm, p[i,t] >= 0)
            @constraint(rm, p[i,t] <= P)
        end
        for i in 1:I
            @constraint(rm, sum(p[i,t] for t=1:T) <= Q)
        end
        # Constraint: cannot exceed inventory limits
        for t in 1:T
            @constraint(rm,
                v1 + sum(p[i,s] for i=1:I, s=1:t) - sum(d[s] for s=1:t) >= Vmin)
            @constraint(rm,
                v1 + sum(p[i,s] for i=1:I, s=1:t) - sum(d[s] for s=1:t) <= Vmax)
        end
    end

    @testset "Affine, manual" begin
        rm = RobustModel(solver=solver)
        # Uncertain parameter: demand at each time stage
        @uncertain(rm, d_nom[t]*(1-θ) <= d[t=1:T] <= d_nom[t]*(1+θ))
        # Decision: how much to produce at each factory at each time
        @variable(rm, p_aff[i=1:I,t=1:T,k=0:t-1])
        p = Any[p_aff[i,t,0] for i in 1:I, t in 1:T]
        for i in 1:I, t in 1:T
            for k in 1:t-1
                p[i,t] += p_aff[i,t,k] * d[k]
            end
        end
        # Objective: total cost of production
        @variable(rm, F); @objective(rm, Min, F)
        @constraint(rm, sum(c[i,t] * p[i,t] for i=1:I, t=1:T) <= F)
        add_constraints(rm, p, d)
        solve(rm)
        @test isapprox(getvalue(F), 44272.82749, atol=TOL)
    end

    @testset "Affine, auto" begin
        rm = RobustModel(solver=solver)
        # Uncertain parameter: demand at each time stage
        @uncertain(rm, d_nom[t]*(1-θ) <= d[t=1:T] <= d_nom[t]*(1+θ))
        # Decision: how much to produce at each factory at each time
        @adaptive(rm, p[i=1:I,t=1:T], policy=Affine, depends_on=d[1:t-1])
        # Objective: total cost of production
        @variable(rm, F); @objective(rm, Min, F)
        @constraint(rm, sum(c[i,t] * p[i,t] for i=1:I, t=1:T) <= F)
        add_constraints(rm, p, d)
        solve(rm)
        @test isapprox(getvalue(F), 44272.82749, atol=TOL)
    end

end  # "with ..."
end  # "Inventory"
