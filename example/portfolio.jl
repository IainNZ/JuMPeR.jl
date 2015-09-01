#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2015: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# example/portfolio.jl
# Solve a robust portfolio optimization problem. Can choose between:
# - Polyhedral or ellipsoidal uncertainty sets
# - Reformulation or cutting planes
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
# Load packages
# Modelling
using JuMP, JuMPeR
# For generating data
using Distributions
# Only needed for Julia 0.3, used for factorizing covariance matrix
using Compat

#-----------------------------------------------------------------------
# Global configuration
# Number of stocks
const NUM_ASSET = 10
# Number of samples of past returns
const NUM_SAMPLES = 1000
# Uncertainty set type (:Polyhedral or :Ellipsoidal)
#const UNCERTAINY_SET = :Polyhedral
const UNCERTAINY_SET = :Ellipsoidal
# Seed the RNG for reproducibilities sake
srand(1234)

#-----------------------------------------------------------------------
# generate_data
# Sample returns of the assets from a multivariate normal distribution.
# Returns a matrix with the samples in rows & each column is an asset.
function generate_data()
    data = zeros(NUM_SAMPLES, NUM_ASSET)

    # Linking factors to induce correlations
    β = [(i-1.0)/NUM_ASSET for i in 1:NUM_ASSET] 

    for i in 1:NUM_SAMPLES
        # Common market factor, μ=3%, σ=5%, truncated at ±3σ
        z = rand(Normal(0.03, 0.05))
        z = max(z, 0.03 - 3*0.05)
        z = min(z, 0.03 + 3*0.05)

        for j in 1:NUM_ASSET
            # Idiosyncratic contribution, μ=0%, σ=5%, truncated at ±3σ
            r = rand(Normal(0.00, 0.05))
            r = max(r, 0.00 - 3*0.05)
            r = min(r, 0.00 + 3*0.05)
            data[i,j] = β[j] * z + r
        end
    end

    return data
end

#-----------------------------------------------------------------------
# solve_portfolio
# Solve the robust portfolio problem given a matrix of past returns.
# The particular problem solved is
#   max   z
#   s.t.  ∑ rᵢ xᵢ ≥ z  ∀ r ∈ R      # Worst-case return
#         ∑ xᵢ    = 1               # Is a portfolio
#         xᵢ ≥ 0
# where R is one of the uncertainty sets:
# R = { (r,z) | r = μ + Σ^½ z,
#        Polyhedral:  ‖z‖₁ ≤ Γ, ‖z‖∞ ≤ 1
#       Ellipsoidal:  ‖z‖₂ ≤ Γ
#     }
function solve_portfolio(past_returns, Γ)
    # Setup the robust optimization model
    m = RobustModel()

    # Each asset is a share of the money to invest...
    @defVar(m, 0 <= x[1:NUM_ASSET] <= 1)
    # ... and we must allocate all the money.
    @addConstraint(m, sum(x) == 1)

    # JuMPeR doesn't support uncertain objectives directly.
    # Uncertain objectives, should put the objective function as
    # a constraint (epigraph form)
    @defVar(m, obj)
    @setObjective(m, Max, obj)

    # The uncertain parameters are the returns for each asset
    @defUnc(m, r[1:NUM_ASSET])

    # The uncertainty set requires the "square root" of the covariance
    μ = vec(mean(past_returns, 1))  # Want a column vector 
    Σ = cov(past_returns)
    L = full(@compat chol(Σ, Val{:L}))  # Σ^½
    # Define auxiliary uncertain parameters to model underlying factors
    @defUnc(m, z[1:NUM_ASSET])
    # Link r with z
    @addConstraint(m, r .== L*z + μ)
    if UNCERTAINY_SET == :Polyhedral
        @addConstraint(m, norm(z,   1) <= Γ)  # ‖z‖₁ ≤ Γ
        @addConstraint(m, norm(z, Inf) <= 1)  # ‖z‖∞ ≤ 1
    else
        @addConstraint(m, norm(z,   2) <= Γ)  # ‖z‖₂ ≤ Γ
    end

    # The objective function: the worst-case return
    @addConstraint(m, obj <= dot(r, x))

    # Solve it
    solve(m)

    # Return the allocation and the worst-case return
    return getValue(x), getValue(obj)
end

#-----------------------------------------------------------------------
# simulate
# Generate some past returns, obtain portfolios for different
# values of Γ, and analyze distribution of returns
function simulate()
    # Generate the simulated data
    past_returns   = generate_data()
    future_returns = generate_data()

    # Print table header
    println("    Γ |    Min |    10% |    20% |   Mean |    Max | Support")
    for Γ in 0:NUM_ASSET
        # Solve for portfolio
        x, obj = solve_portfolio(past_returns, Γ)
        # Evaluate portfolio returns in future
        future_z = future_returns*x
        sort!(future_z)
        min_z    = future_z[1]*100
        ten_z    = future_z[div(NUM_SAMPLES,10)]*100
        twenty_z = future_z[div(NUM_SAMPLES, 5)]*100
        mean_z   = mean(future_z)*100
        max_z    = future_z[end]*100
        support  = sum(x .> 1e-6)  # Number of assets used
        @printf(" %4.1f | %6.2f | %6.2f | %6.2f | %6.2f | %6.2f | %7d\n",
                    Γ, min_z, ten_z, twenty_z, mean_z, max_z, support)
    end
end

simulate()