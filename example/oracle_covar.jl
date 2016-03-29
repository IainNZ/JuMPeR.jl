#############################################################################
# JuMPeR
# Julia for Mathematical Programming - extension for Robust Optimization
# See http://github.com/IainNZ/JuMPeR.jl
#############################################################################
# EXAMPLE
# oracle_covar.jl
# Defines a data-driven oracle that takes a dataset, computes mean and
# covariance information, and uses it to generate cuts.
#
# Given a data matrix D with one row per sample, calculate the mean
# for each column, and the sample covariance matrix S of D. Let A be
# a Cholesky factorization of PSD matrix S = AA', giving the set
#  U = { u | u = A z + μ, ‖z‖ <= Γ }
#
# The cutting plane problem
#  max x'u
# s.t. u ∈ U
# has a closed form solution
#  z = Γ *  A'c / ‖A'c‖
#  u = Az + μ
# with objective value Γ ‖A'c‖ +  c'μ

#############################################################################

using JuMPeR

type CovarOracle <: AbstractOracle
    μ::Vector{Float64}
    S::Matrix{Float64}
    A::Matrix{Float64}
    Γ::Float64
    cut_tol::Float64
end
# Take the raw data and transform it for our needs in the constructor of
# the oracle. We could also do this during setup if we desired
function CovarOracle(raw_data::Matrix{Float64}, Γ::Number)
    S = cov(raw_data)
    return CovarOracle(vec(mean(raw_data, 1)), S, chol(S), float(Γ), 1e-6)
end

# We take no action for individual constraints
function JuMPeR.register_constraint(cv::CovarOracle, rm::Model, ind::Int, prefs)
    nothing
end

# Setup has already been done
function JuMPeR.setup_set(cv::CovarOracle, rm::Model, prefs)
    cv.cut_tol = get(prefs, :cut_tol, 1e-6)
end

# We could add reformulation if wanted to, but we will stick to cutting
# planes here. Reformulations are much more involved as we often need to
# add multiple constraints and variables to the master problem.
function JuMPeR.generate_reform(cv::CovarOracle, master::Model, rm::Model, inds::Vector{Int})
    return 0
end

# Add cutting planes for all constraints where it will change the solution
function JuMPeR.generate_cut(cv::CovarOracle, master::Model, rm::Model, inds::Vector{Int}, active=false)
    # Extract the current master solution
    master_sol = master.colVal

    # Keep track of the constraints that need to be added to the master
    # problem by the main solve loop
    new_cons = {}

    for con_ind in inds
        # Get a reference to the constraint we are looking at
        con = JuMPeR.get_uncertain_constraint(rm, con_ind)

        # Obtain the objective of the cutting plane problem using a built-in
        # utility function in JuMPeR `build_cut_objective`.
        # This function takes a constraint, and returns
        #  - :Min or :Max, depending on the constraint's sign
        #  - The objective coefficients in the cutting problem in the header
        #  - the constant value of the LHS given the current master solution
        cut_sense, c, lhs_const = 
            JuMPeR.build_cut_objective(rm, con, master_sol)

        # Calculate A'c once
        Ac = cv.A'*c

        # Calculate ‖A'c‖ once
        norm_Ac = norm(Ac)

        # Calculate new cut's left-hand-side
        lhs_of_cut = cv.Γ * norm_Ac + dot(c, cv.μ) + lhs_const
        
        # Check violation using JuMPeR utility function `check_cut_status`
        # which, given the constraint in question, the proposed new value
        # of the left-hand-side, and a tolerance, returns either
        # :Violate, :Active, or :Slack
        if JuMPeR.check_cut_status(con, lhs_of_cut, cv.cut_tol) != :Violate
            continue  # No violation, no new cut
        end
        
        # This constraint will violate the new constraint
        # Calculate the value of the uncertainties
        u = cv.A * (cv.Γ * Ac ./ norm_Ac) + cv.μ

        # Use a JuMPeR utility function to take our uncertain constraint
        # and create a deterministic constraint
        new_con = JuMPeR.build_certain_constraint(master, con, u)

        # Add to the list of constraints to return
        push!(new_cons, new_con)
    end
    
    return new_cons
end