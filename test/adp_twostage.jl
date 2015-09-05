#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2015: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# test/adaptive.jl
# Adaptive problems
#-----------------------------------------------------------------------

using JuMP, JuMPeR
using FactCheck

facts("[adaptive] Newsvendor") do

# Generate data
N = 5
srand(N)
Dμ = rand(20:50, N)
Dσ = 0.5 * Dμ

context("Static, manual") do
    m = RobustModel()
    # Budget set
    @defUnc(m, D[1:N])
    @addConstraint(m, norm((D - Dμ)./Dσ,   1) <= isqrt(N))
    @addConstraint(m, norm((D - Dμ)./Dσ, Inf) <= 1)
    # Initial stock bought
    @defVar(m, x >= 0)
    # Adaptive: stock sold
    @defVar(m, S[1:N] >= 0)
    # Max profit
    @defVar(m, z)
    @setObjective(m, Max, z)
    @addConstraint(m, z <= sum(S) - 0.5*x)
    # Can't sell more than demand
    for i in 1:N
        @addConstraint(m, S[i] <= D[i])
    end
    # Can't sell more than bought
    @addConstraint(m, sum(S) <= x)
    solve(m)
    @fact getValue(x) --> roughly(    sum( Dμ-Dσ ), atol=1e-4)
    @fact getValue(z) --> roughly(0.5*sum((Dμ-Dσ)), atol=1e-4)
end

context("Static, auto ") do
    m = RobustModel()
    # Budget set
    @defUnc(m, D[1:N])
    @addConstraint(m, norm((D - Dμ)./Dσ,   1) <= isqrt(N))
    @addConstraint(m, norm((D - Dμ)./Dσ, Inf) <= 1)
    # Initial stock bought
    @defVar(m, x >= 0)
    # Adaptive: stock sold
    @defAdaptVar(m, S[1:N] >= 0, policy=Static, depends_on=D)
    # Max profit
    @defVar(m, z)
    @setObjective(m, Max, z)
    @addConstraint(m, z <= sum(S) - 0.5*x)
    # Can't sell more than demand
    for i in 1:N
        @addConstraint(m, S[i] <= D[i])
    end
    # Can't sell more than bought
    @addConstraint(m, sum(S) <= x)
    solve(m)
    @fact getValue(x) --> roughly(    sum( Dμ-Dσ ), atol=1e-4)
    @fact getValue(z) --> roughly(0.5*sum((Dμ-Dσ)), atol=1e-4)
end




context("Affine, auto ") do
    m = RobustModel()
    # Budget set
    @defUnc(m, D[1:N])
    @addConstraint(m, norm((D - Dμ)./Dσ,   1) <= 1)
    @addConstraint(m, norm((D - Dμ)./Dσ, Inf) <= 1)
    # Initial stock bought
    @defVar(m, x >= 0)
    # Adaptive: stock sold
    @defAdaptVar(m, S[1:N] >= 0, policy=Affine, depends_on=D)
    # Max profit
    @defVar(m, z)
    @setObjective(m, Max, z)
    @addConstraint(m, z <= sum(S) - 0.5*x)
    # Can't sell more than demand
    for i in 1:N
        @addConstraint(m, S[i] <= D[i])
    end
    # Can't sell more than bought
    @addConstraint(m, sum(S) <= x)
    solve(m)
    @show sum(Dμ)
    @show maximum(Dσ)
    @show getValue(x)# --> roughly(    sum( Dμ-Dσ ), atol=1e-4)
    @show getValue(z)# --> roughly(0.5*sum((Dμ-Dσ)), atol=1e-4)
end

end