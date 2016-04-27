#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2016: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# example/inventory.jl
# Implements the inventory management example from
#   A. Ben-Tal, A. Goryashko, E. Guslitzer, A. Nemirovski
#   "Adjustable Robust Solutions of Uncertain Linear Programs"
# Requires a linear optimization solver.
#-----------------------------------------------------------------------

using JuMP, JuMPeR

# Define problem parameters
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

# Setup robust model
invmgmt = RobustModel()

# Uncertain parameter: demand at each time stage lies in a interval
@uncertain(invmgmt, d_nom[t]*(1-θ) <= d[t=1:T] <= d_nom[t]*(1+θ))

# Decision: how much to produce at each factory at each time
# As this decision can be updated as demand is realized, we will use adaptive
# policy - in particular, an affine policy where production at time t is an
# affine function of the demand realized previously.
@adaptive(invmgmt, p[i=1:I,t=1:T], policy=Affine, depends_on=d[1:t-1])

# Objective: minimize total cost of production
@variable(invmgmt, F)  # Overall cost
@objective(invmgmt, Min, F)
@constraint(invmgmt, F >= sum{c[i,t] * p[i,t], i=1:I, t=1:T})

# Constraint: cannot exceed production limits
for i in 1:I, t in 1:T
    @constraint(invmgmt, p[i,t] >= 0)
    @constraint(invmgmt, p[i,t] <= P)
end
for i in 1:I
    @constraint(invmgmt, sum{p[i,t], t=1:T} <= Q)
end

# Constraint: cannot exceed inventory limits
for t in 1:T
    @constraint(invmgmt,
        v1 + sum{p[i,s], i=1:I, s=1:t} - sum{d[s],s=1:t} >= Vmin)
    @constraint(invmgmt,
        v1 + sum{p[i,s], i=1:I, s=1:t} - sum{d[s],s=1:t} <= Vmax)
end

# Solve
status = solve(invmgmt)

println(getobjectivevalue(invmgmt))
