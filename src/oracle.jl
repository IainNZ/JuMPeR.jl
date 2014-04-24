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

# build_certain_constraint                              [not exported, pure]
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

#############################################################################
# Default included oracles
include("oracle_poly.jl")         # PolyhedralOracle
include("oracle_bertsim.jl")      # BertSimOracle