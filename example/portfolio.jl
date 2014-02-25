#############################################################################
# JuMPeR
# Julia for Mathematical Programming - extension for Robust Optimization
# See http://github.com/IainNZ/JuMPeR.jl
#############################################################################
# portfolio.jl
# Recreates the results from the paper "Data-Driven Robust Optimization" by
# Bertsimas, Gupta, and Kallus, 2014.
#############################################################################

using JuMPeR
import JuMPeR: registerConstraint, setup, generateCut, generateReform
using Distributions

#############################################################################
# SetMOracle
# In text                   U^M
# Structural assumptions    None
# Hypothesis test           Marginal samples
# Note that the oracle is setup very precisely for this problem, in that
# trying to use it on multiple constraints isn't going to work. It wouldn't
# be a big deal to change that, although presumably you'd have a different
# data vector for each constraint so you'd probably want a different oracle
# for each on anyway.
#############################################################################
type SetMOracle <: AbstractOracle
    con
    epsilon::Float64
    delta::Float64
    data
    s::Int
end
SetMOracle(eps, delta, data) = SetMOracle(nothing, eps, delta, data, 0)

# We only do reformulations - could do cuts if wanted to support x being
# nonnegative (because would have something to do - put at upper or lower
# bound)
function registerConstraint(w::SetMOracle, con, ind::Int, prefs)
    w.con = con
    return [:Cut => false, :Reform => true, :Sample => false]
end

# This could be done in constructor too, since we don't use any properties
# of the model. But we'll leave it here for now.
function setup(w::SetMOracle, rm::Model)
    # Calculate s from Eq (28)
    N = size(w.data, 1)  # Number of samples
    d = size(w.data, 2)  # Dimension of vector
    binom_dist = Distributions.Binomial(N, 1 - w.epsilon/d)
    w.s = int(quantile(binom_dist, 1 - w.delta / (2*d)) + 1)
            # Need +1 to be consistent with indexing in paper - could refactor
    println("SetMOracle.s = $(w.s) for N = $N")

    # Sort marginals in ascending order (makes a copy... efficient
    # sort-in-place for data stored in matrix like this not possible in Julia
    # v0.2 but will be maybe as soon as v0.3 using "views")
    for dim = 1:d
        w.data[1:end, dim] = sort(w.data[1:end, dim])
    end
end

# Not used - see registerConstraint
generateCut(w::SetMOracle, rm::Model, ind::Int, m::Model) = return 0

# Where we apply the set in Eq (31)
# TODO: Make less brittle, assumes that sense of constraint is exactly as it
#       is in this example code, so sets everything to lower bound.
function generateReform(w::SetMOracle, rm::Model, ind::Int, m::Model)
    N = size(w.data, 1)  # Number of samples
    d = size(w.data, 2)  # Dimension of vector
    s = w.s

    if N - s + 1 >= d
        println("Condition from paper not met for probabilistic guarantee! N=$N, s=$s")
    end

    # We need need to build the new constraints
    # We make a new "affine expression" for the left-hand-side that consists
    # only of the certain part of the "original" left-hand-side
    orig_lhs = w.con.terms
    new_lhs = AffExpr(orig_lhs.vars,
                      [orig_lhs.coeffs[i].constant::Float64 
                                for i in 1:length(orig_lhs.vars)],
                      orig_lhs.constant.constant)

    for var_ind = 1:length(orig_lhs.vars)
        num_uncs = length(orig_lhs.coeffs[var_ind].uncs)
        if num_uncs == 0
            # Not one of the coefficients we are interested in
            # That is, the variable is z
            continue
        elseif num_uncs > 1
            # What is this?!
            error("Only designed for one coefficient per x")
        end
        unc = orig_lhs.coeffs[var_ind].uncs[1].unc  # Use as asset index
        new_lhs.coeffs[var_ind] += w.data[(N-s+1)+1, unc]  # Correct for 0-base
    end

    # Generate the "reformulation"
    @addConstraint(m, new_lhs >= w.con.lb)

    return true
end

#############################################################################
# End of SetMOracle
#############################################################################

#############################################################################
# Simulate the data
#############################################################################
function generate_data(N)
    print("Generating data...")
    d    = 10                       # Number of assets
    beta = [(i-1.)/9. for i = 1:10] # Linking factors
    data = zeros(N, d)

    for sample_ind = 1:N
        # Common market factor, mean 3%, sd 5%, truncate at +- 3 sd
        z = rand(Normal(0.03, 0.05))
        z = max(z, 0.03 - 3*0.05)
        z = min(z, 0.03 + 3*0.05)

        for asset_ind = 1:d
            # Idiosyncratic contribution, mean 0%, sd 5%, truncated at +- 3 sd
            asset = rand(Normal(0.00, 0.05))
            asset = max(asset, 0.00 - 3*0.05)
            asset = min(asset, 0.00 + 3*0.05)
            data[sample_ind, asset_ind] = beta[asset_ind] * z + asset
        end
    end

    println(" done")
    for asset_ind = 1:d
        println("Min value for asset $asset_ind is $(minimum(data[1:end,asset_ind]))")
    end
    return data
end

#############################################################################
# Solve the robust portfolio problem
#############################################################################
function solve_portfolio(past_returns)

    # Uncertainty set parameters
    epsilon = 0.10
    delta   = 0.20

    # Create oracle
    oracle = SetMOracle(epsilon, delta, past_returns)

    # Setup the robust optimization model
    m = RobustModel()
    @defVar(m, z)
    @defVar(m, x[1:10] >= 0)
    @defUnc(m, r[1:10])

    @setObjective(m, Max, z)
    addConstraint(m, sum([      x[i] for i=1:10 ])     == 1)
    addConstraint(m, sum([ r[i]*x[i] for i=1:10 ]) - z >= 0, oracle)

    # Solve it, report statistics on number of cutting planes etc.
    println("Solving model...")
    solveRobust(m, report=true)

    # Solution
    println(getValue(x))
    println(getValue(z))
end

past_returns = generate_data(1200)
solve_portfolio(past_returns)

past_returns = generate_data(1200)
println("\nRunning again to see how long it 'really' takes")
solve_portfolio(past_returns)