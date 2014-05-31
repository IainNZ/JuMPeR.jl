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
# passes any preferences provided via the solve command. Returns a dictionary
# where the keys are the symbols :Cut, :Reform, and :Sample and the values
# are true or false, where true indicates the oracle has selected these
# operations for this constraint. (not currently used)
registerConstraint(w::AbstractOracle, con, ind::Int, prefs) = error("Not implemented!")

# setup
# Gives oracle time to do any setup it needs to do. Called after all
# constraints have been registered. Examples of work that could be done here
# include transforming the uncertainty set and generating a cutting plane
# model. May be called multiple times - this should be handled by the oracle
setup(w::AbstractOracle, rm::Model) = error("Not implemented")

# generateCut
# Called in the main loop every iteration. m is the actual current model, aka
# the master model, that will have the current solution and to which 
# constraints should be added. Returns number of cuts added.
# cb is nothing iff problem contains no integer variables, will be the
# callback handle otherwise. If provided, should be used to add lazy 
# constraints instead of normal constraints.
generateCut(w::AbstractOracle, rm::Model, ind::Int, m::Model, cb) = error("Not implemented")

# generateReform
# Called before the main loop, adds anything it wants to the model
generateReform(w::AbstractOracle, rm::Model, ind::Int, m::Model) = error("Not implemented")

export AbstractOracle
export registerConstraint, setup, generateCut, generateReform
export setDefaultOracle

function setDefaultOracle!(rm, w::AbstractOracle)
    getRobust(rm).defaultOracle = w
end

#############################################################################
# Helper functions that can be shared by all oracles

# build_certain_constraint
# Takes an uncertain constraint (unc_con) that we are making certain by
# replacing uncertainties with actual numbers (unc_val) - the values of all
# uncertainties in the robust optimization problem
function build_certain_constraint(  unc_con::UncConstraint, 
                                    unc_val::Vector{Float64} )
    unc_lhs = unc_con.terms   
    num_var = length(unc_lhs.vars)
    new_lhs = AffExpr(unc_lhs.vars,
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

function build_certain_constraint(  unc_con::UncConstraint, 
                                    unc_val_dict::JuMPDict{Float64} )
    length(unc_val_dict.indexsets) != 1 &&
        error("only JuMPDict with one dimension allowed. Manually convert to Vector{Float64}")
    return build_certain_constraint(unc_con, vec(unc_val_dict.innerArray))
end

# build_cut_objective
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
function build_cut_objective(   unc_con::UncConstraint,
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