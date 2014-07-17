#############################################################################
# JuMPeR
# Julia for Mathematical Programming - extension for Robust Optimization
# See http://github.com/IainNZ/JuMPeR.jl
#############################################################################
# portfolio.jl
# Use the CovarOracle defined in oracle_covar.jl to solve an asset
# allocation problem. Solves for multiple values of Γ and shows the
# distribution of returns of the optimal portfolio for each value.
#############################################################################

using JuMPeR
using Distributions
using Gurobi
include("oracle_covar.jl")

const NUM_ASSET = 10

# Seed the RNG for reproducibilities sake
srand(10)

#############################################################################
# generate_data
# Sample returns of the assets. The parameters of the assets are fixed, we
# simply make num_samples draws from the joint distribution.
# Returns matrix with the samples in rows, each column is an assets
function generate_data(num_samples::Int)
    data = zeros(num_samples, NUM_ASSET)

    # Linking factors
    beta = [(i-1.0)/NUM_ASSET for i = 1:NUM_ASSET] 

    for sample_ind = 1:num_samples
        # Common market factor, mean 3%, sd 5%, truncate at +- 3 sd
        z = rand(Normal(0.03, 0.05))
        z = max(z, 0.03 - 3*0.05)
        z = min(z, 0.03 + 3*0.05)

        for asset_ind = 1:NUM_ASSET
            # Idiosyncratic contribution, mean 0%, sd 5%, truncated at +- 3 sd
            asset = rand(Normal(0.00, 0.05))
            asset = max(asset, 0.00 - 3*0.05)
            asset = min(asset, 0.00 + 3*0.05)
            data[sample_ind, asset_ind] = beta[asset_ind] * z + asset
        end
    end

    return data
end


#############################################################################
# solve_portfolio
# Solve the robust portfolio problem given a matrix of past returns and
# a degree of conservatism, Γ. The particular problem solved is
#   max   z
#   s.t.  ∑ xi    = 1
#         ∑ ri*xi ≥ z  ∀ r ∈ R
#         xi ≥ 0
#   where R is CovarOracle's set, an ellipse centered on the mean return
#         and tilted using the covariance of the returns
function solve_portfolio(past_returns, Γ)
    # Setup the robust optimization model
    m = RobustModel(solver=GurobiSolver(OutputFlag=0))

    # Create the CovarOracle
    setDefaultOracle!(m, CovarOracle(past_returns, Γ))

    # Variables
    @defVar(m, obj)  # Put objective as constraint (epigraph form)
    @defVar(m, 0 <= x[1:NUM_ASSET] <= 1)  # Fractional allocation

    # Uncertainties
    @defUnc(m, r[1:NUM_ASSET])  # The returns

    @setObjective(m, Max, obj)

    # Portfolio constraint
    @addConstraint(m, sum(x) == 1)

    # The objective constraint - uncertain
    @addConstraint(m, obj - dot(r, x) <= 0)

    # Solve it, report statistics on number of cutting planes etc.
    #println("Solving model for Γ=$Γ")
    #solveRobust(m, report=true)
    #println("Objective value: ", getValue(obj))
    solveRobust(m)

    return getValue(x), getValue(obj)
end


#############################################################################
# simulate
# Generate some past returns, then obtain portfolios for different values
# of Γ to see the effect of robustness on the distribution of returns
function simulate(num_past, num_future)
    # Generate the simulated data
    past_returns   = generate_data(num_past)
    future_returns = generate_data(num_future)

    # Print table header
    println("   Γ|   Min|   10%|   20%|  Mean|   Max")
    for Γ in [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
        # Solve for portfolio
        x, obj = solve_portfolio(past_returns, Γ)
        # Evaluate portfolio returns in future
        future_z = future_returns*x[:]
        sort!(future_z)
        min_z    = future_z[1]*100
        ten_z    = future_z[int(num_future*0.1)]*100
        twenty_z = future_z[int(num_future*0.2)]*100
        mean_z   = mean(future_z)*100
        max_z    = future_z[end]*100
        @printf(" %3.1f| %5.1f| %5.1f| %5.1f| %5.1f| %5.1f\n",
                    Γ, min_z, ten_z, twenty_z, mean_z, max_z)
    end
end

simulate(1000, 1000)