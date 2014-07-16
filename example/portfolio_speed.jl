#############################################################################
# JuMPeR
# Julia for Mathematical Programming - extension for Robust Optimization
# See http://github.com/IainNZ/JuMPeR.jl
#############################################################################
# portfolio_speed.jl
# A counter-part to the portfolio.jl example. Compares the speed difference
# between a dedicated oracle (CovarOracle) and the general oracle that uses
# Gurobi to solve the cutting plane problem
#############################################################################

using JuMPeR
using Distributions
using Gurobi
include("oracle_covar.jl")

const NUM_ASSET = 20

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
function solve_portfolio(past_returns, Γ, oracle_choice, report_time)
    # Setup the robust optimization model
    m = RobustModel(solver=GurobiSolver(OutputFlag=0))

    # Variables
    @defVar(m, obj)  # Put objective as constraint (epigraph form)
    @defVar(m, 0 <= x[1:NUM_ASSET] <= 1)  # Fractional allocation

    # Uncertainties
    @defUnc(m, r[1:NUM_ASSET])  # The returns

    @setObjective(m, Max, obj)

    # Portfolio constraint
    addConstraint(m, sum(x) == 1)

    # The objective constraint - uncertain
    addConstraint(m, obj - dot(r, x) <= 0)

    if oracle_choice == :CovarOracle
        # Create the CovarOracle
        setDefaultOracle!(m, CovarOracle(past_returns, Γ))
    elseif oracle_choice == :GeneralOracleCut ||
           oracle_choice == :GeneralOracleReform
        # Setup the uncertainty set manually
        μ = vec(mean(past_returns, 1))
        S = cov(past_returns)
        A = chol(S)
        @defUnc(m, z[1:NUM_ASSET])
        for i = 1:NUM_ASSET
            addConstraint(m, r[i] == sum([A[i,j]*z[j] for j=1:NUM_ASSET]) + μ[i])
        end
        addEllipseConstraint(m, [z[i] for i=1:NUM_ASSET], Γ)
    end

    # Solve it
    solveRobust(m,  prefer_cuts = (oracle_choice == :GeneralOracleCut),
                    report      = report_time,
                    cut_tol     = 1e-5)

    return getValue(x), getValue(obj)
end


#############################################################################
# simulate
# Generate some past returns, then solve the portfolio problem using first
# a dedicated oracle, then the general oracle included in JuMPeR (with a
# warm-up run first to remove transient effects). The results show that
# the reformulation using the GeneralOracle is best in this case, followed
# by the CovarOracle's cutting planes and finally the GeneralOracle's
# cutting planes.
function simulate(num_past, num_future)
    # Generate the simulated data
    past_returns   = generate_data(num_past)
    
    # Set the conservatism level
    Γ = 3.0

    # Warmup runs
    x, obj = solve_portfolio(past_returns, Γ, :CovarOracle, false)
    x, obj = solve_portfolio(past_returns, Γ, :GeneralOracleCut, false)
    x, obj = solve_portfolio(past_returns, Γ, :GeneralOracleReform, false)

    # Solve with CovarOracle
    println("== CovarOracle          ====================")
    x, obj = solve_portfolio(past_returns, Γ, :CovarOracle, true)
    println("Objective value: ", obj)

    # Solve with GeneralOracle
    println("== GeneralOracle-Cuts   ====================")
    x, obj = solve_portfolio(past_returns, Γ, :GeneralOracleCut, true)
    println("Objective value: ", obj)
    println("== GeneralOracle-Reform ====================")
    x, obj = solve_portfolio(past_returns, Γ, :GeneralOracleReform, true)
    println("Objective value: ", obj)
end

simulate(1000, 1000)