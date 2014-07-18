#############################################################################
# JuMPeR
# Julia for Mathematical Programming - extension for Robust Optimization
# See http://github.com/IainNZ/JuMPeR.jl
#############################################################################
# Oracles
# Oracles are "robustifying operators". They take a list of constraints,
# some of which maybe null or partially defined, and are reponsible for
# "robustifying" them by either reformulating them, providing a cutting plane
# algorithm, and/or a sampling procedure. They have the full statement of
# the uncertain problem available to them, in particular the uncertainty set
# and uncertain bounds.
#############################################################################

#############################################################################
# AbstractOracle
# All oracles implement the interface defined by AbstractOracle
abstract AbstractOracle

# registerConstraint
# Notifies the oracle that it is responsible for this constraint, and 
# passes any preferences provided via the solveRobust() command.
registerConstraint(ab::AbstractOracle, rm::Model, ind::Int, prefs) = 
    error("$(typeof(ab)) hasn't implemented registerConstraint")

# setup
# Gives oracle time to do any setup it needs to do. Called after all
# constraints have been registered. Examples of work that could be done here
# include transforming the uncertainty set and generating a cutting plane
# model. Will NOT be called multiple times.
setup(ab::AbstractOracle, rm::Model, prefs) = 
    error("$(typeof(ab)) hasn't implemented setup")

# generateReform
# Called before the main loop, adds anything it wants to the model. Returns
# number of constraints reformulated.
generateReform(ab::AbstractOracle, master::Model, rm::Model, inds::Vector{Int}) =
    error("$(typeof(ab)) hasn't implemented generateReform")

# generateCut
# Called in the main loop every iteration/every time an integer solution is
# found. Returns a vector of constraints which are added to the problem by
# the main solve loop.
# The optional "active" argument will be called if the user wants to know
# the active scenarios at optimality.
generateCut(ab::AbstractOracle, master::Model, rm::Model, inds::Vector{Int}, active=false) =
    error("$(typeof(ab)) hasn't implemented generateCut")


export AbstractOracle
export registerConstraint, setup, generateCut, generateReform
export setDefaultOracle

function setDefaultOracle!(rm, w::AbstractOracle)
    getRobust(rm).defaultOracle = w
end

#############################################################################
# Utility functions that can be shared by all oracles

# get_uncertain_constraint
# Given a RobustModel and a constraint index, returns the uncertain
# constraint for that index
get_uncertain_constraint(rm::Model, ind::Int) =
    getRobust(rm).uncertainconstr[ind]

# build_certain_constraint
# Takes an uncertain constraint (unc_con) that we are making certain by
# replacing uncertainties with actual numbers (unc_val) - the values of all
# uncertainties in the robust optimization problem
function build_certain_constraint(  master::Model,
                                    unc_con::UncConstraint, 
                                    unc_val::Vector{Float64} )
    unc_lhs = unc_con.terms
    num_var = length(unc_lhs.vars)
    new_lhs = AffExpr([Variable(master,v.col) for v in unc_lhs.vars],
                      [unc_lhs.coeffs[i].constant for i in 1:num_var],
                      0.0)
    
    # Variable part
    for var_ind = 1:num_var
        coeff = unc_lhs.coeffs[var_ind]
        for unc_ind = 1:length(coeff.vars)
            new_lhs.coeffs[var_ind] += unc_val[coeff.vars[unc_ind].unc] *
                                               coeff.coeffs[unc_ind]
        end
    end
    # Non variable part
        coeff = unc_lhs.constant
        for unc_ind = 1:length(coeff.vars)
            new_lhs.constant += unc_val[coeff.vars[unc_ind].unc] *
                                        coeff.coeffs[unc_ind]
        end

    return sense(unc_con) == :(<=) ? new_lhs <= unc_con.ub :
                                     new_lhs >= unc_con.lb
end

function build_certain_constraint(  master::Model,
                                    unc_con::UncConstraint, 
                                    unc_val_dict::JuMPDict{Float64} )
    length(unc_val_dict.indexsets) != 1 &&
        error("only JuMPDict with one dimension allowed. Manually convert to Vector{Float64}")
    return build_certain_constraint(master, unc_con, vec(unc_val_dict.innerArray))
end

# build_cut_objective
# Takes an uncertain constraint (unc_con) and a master solution (x_val)
# and returns the coefficients for each uncertainty, as
# well as the objective sense and the constant term of the objective
# of the cutting plane problem that arises from variables that do not have
# uncertain coefficients.
# For example, for input:
#     unc_con = (3*u[1] + 2.0) * x[1] + (  u[2] - 1.0) * x[2] +
#               (u[1] +  u[3]) * x[3] + (u[3] +2*u[4]) * x[4] <= 5.0 + u[5]
#     x_val   = [2.0, 3.0, 4.0, 5.0]
# returns
#     :Max, [], 1.0
# Constraint with build_cut_objective_sparse, which returns only
# the coefficients that appear in the constraint
function build_cut_objective(   rm::Model,
                                unc_con::UncConstraint,
                                x_val::Vector{Float64})
    unc_coeffs = zeros(getNumUncs(rm))
    unc_lhs    = unc_con.terms
    lhs_const  = 0.0

    # Uncertains attached to variables
    for var_ind = 1:length(unc_lhs.vars)
        uaff = unc_lhs.coeffs[var_ind]
        col  = unc_lhs.vars[var_ind].col
        for unc_ind = 1:length(uaff.vars)
            unc = uaff.vars[unc_ind].unc
            unc_coeffs[unc] += uaff.coeffs[unc_ind] * x_val[col]
        end
        lhs_const += uaff.constant * x_val[col]
    end
    # Uncertains not attached to variables
        uaff = unc_lhs.constant
        for unc_ind = 1:length(uaff.vars)
            unc = uaff.vars[unc_ind].unc
            unc_coeffs[unc] += uaff.coeffs[unc_ind]
        end

    return (sense(unc_con) == :(<=) ? :Max : :Min), 
                unc_coeffs, lhs_const
end


# build_cut_objective_sparse
# Takes an uncertain constraint (unc_con) and a master solution (x_val)
# and returns the coefficients for each uncertain in the cutting plane, as
# well as the objective sense and the constant term of the objective
# of the cutting plane problem that arises from variables that do not have
# uncertain coefficients.
# For example, for input:
#     unc_con = (3*u[1] + 2.0) * x[1] + (  u[2] - 1.0) * x[2] +
#               (u[1] +  u[3]) * x[3] + (u[3] +2*u[4]) * x[4] <= 5.0 + u[5]
#     x_val   = [2.0, 3.0, 4.0, 5.0]
# returns
#     :Max, [(5,-1.0),(4,10.0),(2,3.0),(3,9.0),(1,10.0)], 1.0
# Note that exact order of coefficients is non-deterministic due to the
# use of a dictionary internally.
function build_cut_objective_sparse(   unc_con::UncConstraint,
                                x_val::Vector{Float64})
    unc_coeffs = Dict{Int,Float64}()
    unc_lhs = unc_con.terms   
    num_var = length(unc_lhs.vars)
    lhs_constant = 0.0

    # Uncertains attached to variables
    for var_ind = 1:num_var
        uaff = unc_lhs.coeffs[var_ind]
        col  = unc_lhs.vars[var_ind].col
        for unc_ind = 1:length(uaff.vars)
            unc = uaff.vars[unc_ind].unc
            if !haskey(unc_coeffs, unc)
                unc_coeffs[unc]  = uaff.coeffs[unc_ind] * x_val[col]
            else
                unc_coeffs[unc] += uaff.coeffs[unc_ind] * x_val[col]
            end
        end
        lhs_constant += uaff.constant * x_val[col]
    end
    # Uncertains not attached to variables
        uaff = unc_lhs.constant
        for unc_ind = 1:length(uaff.vars)
            unc = uaff.vars[unc_ind].unc
            if !haskey(unc_coeffs, unc)
                unc_coeffs[unc]  = uaff.coeffs[unc_ind]
            else
                unc_coeffs[unc] += uaff.coeffs[unc_ind]
            end
        end

    return (sense(unc_con) == :(<=) ? :Max : :Min), 
                collect(unc_coeffs), lhs_constant
end

# check_cut_status
# A simple helper function to check whether a left-hand-side will violate
# an inequality constraint, is binding, or is loose
# Returns  :Slack   -  more than tol gap between lhs and rhs
#          :Active  -  absolute gap between lhs and rhs < gap
#          :Violate -  constraint violated by more than tol
function check_cut_status(con, lhs_value, tol)
    abs(lhs_value - rhs(con)) <= tol && return :Active
    if sense(con) == :<=
        lhs_value > con.ub + tol && return :Violate    
    else
        @assert sense(con) == :>=
        lhs_value < con.lb - tol && return :Violate
    end
    return :Slack
end

#############################################################################
# Default included oracles
include("oracle_gen.jl")            # GeneralOracle
include("oracle_bertsim.jl")        # BertSimOracle