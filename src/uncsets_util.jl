#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2016: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# src/uncsets_util.jl
# Various utility functions to aid in writing uncertainty sets.
# Defines:
#   check_cut_status
# Included by src/uncsets.jl
#-----------------------------------------------------------------------


"""
    check_cut_status(UncConstraint, lhs_value, tol)

Given an uncertain constraint UncConstraint, and a proposed value for the
left-hand-side (typically generated while attempting to find a new cut),
determine the status of the new constraint. Returns one of (for <=):
  * `:Slack` - `rhs - lhs > tol`
  * `:Active` - `|rhs - lhs| <= tol`
  * `:Violate` - `lhs - rhs > tol`
"""
function check_cut_status(con, lhs_value::Float64, tol::Float64)
    gap = rhs(con) - lhs_value
    abs(gap) <= tol && return :Active
    if sense(con) == :(<=)
        # lhs <= rhs  -> lhs - rhs <= tol
        lhs_value - rhs(con) > tol && return :Violate
    else
        # lhs >= rhs  -> rhs - lhs <= tol
        rhs(con) - lhs_value > tol && return :Violate
    end
    return :Slack
end


"""
    build_certain_constraint(UncConstraint, unc_val::Vector{Float64})

Build a deterministic constraint from an uncertain constraint given values
for the uncertain parameters.
"""
function build_certain_constraint(unc_con::UncConstraint, unc_val::Vector{Float64})
    new_lhs = AffExpr()
    function unc_expr_to_coeff(unc_expr)
        new_coeff = unc_expr.constant
        for (unc_coeff, unc) in linearterms(unc_expr)
            new_coeff += unc_coeff * unc_val[unc.id]
        end
        return new_coeff
    end
    # Uncertain expression--variable terms
    for (unc_expr, var) in linearterms(unc_con.terms)
        push!(new_lhs, unc_expr_to_coeff(unc_expr), var)
    end
    # Standalone uncertain expression/constant term
    new_lhs.constant = unc_expr_to_coeff(unc_con.terms.constant)
    # Build the constraint
    if sense(unc_con) == :(<=)
        return @LinearConstraint(new_lhs <= rhs(unc_con))
    else
        return @LinearConstraint(new_lhs >= rhs(unc_con))
    end
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
    rme = get_robust(rm)
    unc_coeffs = zeros(rme.num_uncs)
    unc_lhs    = unc_con.terms
    lhs_const  = 0.0

    # Uncertains attached to variables
    for var_ind = 1:length(unc_lhs.vars)
        uaff = unc_lhs.coeffs[var_ind]
        col  = unc_lhs.vars[var_ind].col
        for unc_ind = 1:length(uaff.vars)
            unc = uaff.vars[unc_ind].id
            unc_coeffs[unc] += uaff.coeffs[unc_ind] * x_val[col]
        end
        lhs_const += uaff.constant * x_val[col]
    end
    # Uncertains not attached to variables
        uaff = unc_lhs.constant
        for unc_ind = 1:length(uaff.vars)
            unc = uaff.vars[unc_ind].id
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
            unc = uaff.vars[unc_ind].id
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
            unc = uaff.vars[unc_ind].id
            if !haskey(unc_coeffs, unc)
                unc_coeffs[unc]  = uaff.coeffs[unc_ind]
            else
                unc_coeffs[unc] += uaff.coeffs[unc_ind]
            end
        end

    return (sense(unc_con) == :(<=) ? :Max : :Min),
                collect(unc_coeffs), lhs_constant
end




"""
    aff_to_uaff(UncExpr, Vector{Variable})

Utility function for uncertainty sets. Given a `UncExpr` and a list of
variables, create an `AffExpr` such that `Uncertain(i)` maps to `Variable(i)`
in the new expression, where `i` is the index.

For example:
    aff_to_uaff(5u[1] + 2u[3], [x[1],x[2],x[3]])  ⟶  5x[1] + 2x[3]
"""
uaff_to_aff(uaff::UncExpr, x::Vector{Variable}) =
    AffExpr(Variable[x[up.id] for up in uaff.vars],
            copy(uaff.coeffs), uaff.constant)
