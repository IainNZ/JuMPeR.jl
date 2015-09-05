#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2015: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# example/inventory.jl
# Implements the inventory management example from
#   A. Ben-Tal, A. Goryashko, E. Guslitzer, A. Nemirovski
#   "Adjustable Robust Solutions of Uncertain Linear Programs"
#-----------------------------------------------------------------------

using JuMP, JuMPeR
using Distributions

#-----------------------------------------------------------------------
# Data
I = 3       # Number of factories
T = 24      # Number of time periods
# Nominal demand
d_nom = 1000*[1 + 0.5*sin(π*(t-1)/12)  for t in 1:T]
θ = 0.20    # Uncertainty level
# Production costs
α = [1.0, 1.5, 2.0]
c = [α[i] * (1 + 0.5*sin(π*(t-1)/12))  for i in 1:I, t in 1:T]
P = 567     # Maximimum production per period
Q = 13600   # Maximumum production over all
Vmin = 500    # Minimum inventory at warehouse
Vmax = 2000   # Maximum inventory at warehouse
v1 = Vmin   # Initial inventory (not provided in paper)


#-----------------------------------------------------------------------
# Setup robust model
invmgmt = RobustModel()

# Uncertain parameter: demand at each time stage
@defUnc(invmgmt, d_nom[t]*(1-θ) <= d[t=1:T] <= d_nom[t]*(1+θ))
# Decision is how much to produce at each factory at each time
# We are using an affine decision rule: the decision is a function of
# the demand already realized, i.e.
# p_i(t) = π_i,t^0  +  ∑_r  π_i,t^r  d_r
@defAdaptVar(invmgmt, p[i=1:I,t=1:T], policy=Affine, depends_on=d[1:t-1])

# Objective: total cost of production
@defVar(invmgmt, F)  # Overall cost
@setObjective(invmgmt, Min, F)
@addConstraint(invmgmt, sum{c[i,t] * p[i,t], i=1:I, t=1:T} <= F)
# Constraint: cannot exceed production limits
for i in 1:I, t in 1:T
    @addConstraint(invmgmt, p[i,t] >= 0)
    @addConstraint(invmgmt, p[i,t] <= P)
end
for i in 1:I
    @addConstraint(invmgmt, sum{p[i,t], t=1:T} <= Q)
end
# Constraint: cannot exceed inventory limits
for t in 1:T
    @addConstraint(invmgmt,
        v1 + sum{p[i,s], i=1:I, s=1:t} - sum{d[s],s=1:t} >= Vmin)
    @addConstraint(invmgmt,
        v1 + sum{p[i,s], i=1:I, s=1:t} - sum{d[s],s=1:t} <= Vmax)
end

# Solve
status = solve(invmgmt)
@show status

#-----------------------------------------------------------------------
# Simulation
# Disabled until I have a way to get the policy back out
#=
# Store the affine policy
π_val = getValue(π)
# Compare affine with deterministic
affine_objs = Float64[]
determ_objs = Float64[]
for sim in 1:1000
    # Draw a sample from the uncertainty set
    d_hat = d_nom .* [rand(Uniform(1-θ,1+θ)) for t in 1:T]
    # Evaluate cost with robust policy
    phat = zeros(I,T)
    for i in 1:I, t in 1:T
        phat[i,t] = π_val[i,t,0]
        if t > 1 
            phat[i,t] += sum([ π_val[i,t,r] * d_hat[r]  for r in 1:t-1 ])
        end
    end
    push!(affine_objs, sum(c .* phat))
    # Get objective for deterministic problem
    det = Model()
    @defVar(det, 0 <= p[i=1:I,t=1:T] <= P)
    @defVar(det, Vmin <= v[t=2:T] <= Vmax)
    @setObjective(det, Min, sum{c[i,t] * p[i,t], i=1:I, t=1:T})
    for i in 1:I
        @addConstraint(det, sum{p[i,t], t=1:T} <= Q)
    end
    @addConstraint(det, v[2] == v1 + sum{p[i,1],i=1:I} - d_hat[1])
    for t in 2:T-1
        @addConstraint(det, v[t+1] == v[t] + sum{p[i,t],i=1:I} - d_hat[t])
    end
    @addConstraint(det, v[T] + sum{p[i,T],i=1:I} - d_hat[T] >= Vmin)
    @addConstraint(det, v[T] + sum{p[i,T],i=1:I} - d_hat[T] <= Vmax)
    solve(det)
    push!(determ_objs, getObjectiveValue(det))
end

@show mean(affine_objs)
@show mean(determ_objs)
@show std(affine_objs)
@show std(determ_objs)
cvar(data,α=0.20) = mean(sort(data)[1:floor(Int,length(data)*α)])
@show cvar(affine_objs)
@show cvar(determ_objs)
=#
#=
using UnicodePlots
rng = linspace(minimum(determ_objs), maximum(affine_objs), 10)
labs, vals = hist(affine_objs, rng)
barplot(collect(labs), vals, title="Affine")
labs, vals = hist(determ_objs, rng)
barplot(collect(labs), vals, title="Optimal")
=#