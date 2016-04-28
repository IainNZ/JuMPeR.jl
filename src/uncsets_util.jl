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
    gap = JuMP.rhs(con) - lhs_value
    abs(gap) <= tol && return :Active
    if JuMP.sense(con) == :(<=)
        # lhs <= rhs  -> lhs - rhs <= tol
        lhs_value - JuMP.rhs(con) > tol && return :Violate
    else
        # lhs >= rhs  -> rhs - lhs <= tol
        JuMP.rhs(con) - lhs_value > tol && return :Violate
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
    if JuMP.sense(unc_con) == :(<=)
        return @LinearConstraint(new_lhs <= JuMP.rhs(unc_con))
    else
        return @LinearConstraint(new_lhs >= JuMP.rhs(unc_con))
    end
end


"""
    build_cut_objective_dense(RobustModel, unc_con)

Takes an uncertain constraint (unc_con) belonging to RobustModel. Returns the
objective function coefficients for the uncertain parameters for a cutting
plane problem for this constraint, using the current RobustModel solution.
Also returns the sense of that cutting plane problem, based on the sense of
the constraint, and the constant term corresponding to the deterministic part
of the constraint.
"""
function build_cut_objective_dense(rm::Model, unc_con::UncConstraint)
    rmext      = get_robust(rm)
    unc_coeffs = zeros(rmext.num_uncs)
    lhs_const  = 0.0
    con_sense  = (JuMP.sense(unc_con) == :(<=)) ? :Max  :  :Min
    # Uncertain expression--variable terms
    for (unc_expr, var) in linearterms(unc_con.terms)
        var_val = getvalue(var)
        for (unc_coeff, unc) in linearterms(unc_expr)
            unc_coeffs[unc.id] += unc_coeff * var_val
        end
        lhs_const += unc_expr.constant * var_val
    end
    # Standalone uncertain expression/constant term
        for (unc_coeff, unc) in linearterms(unc_con.terms.constant)
            unc_coeffs[unc.id] += unc_coeff
        end
        lhs_const += unc_con.terms.constant.constant
    return con_sense, unc_coeffs, lhs_const
end


"""
    build_cut_objective_sparse(RobustModel, unc_con)

Takes an uncertain constraint (unc_con) belonging to RobustModel. Returns the
objective function coefficients for the uncertain parameters for a cutting
plane problem for this constraint, using the current RobustModel solution.
Unlike `build_cut_objective_dense`, the return is sparse: the coefficients
are returned as a vector of tuples (uncertain.id, coeff).
Also returns the sense of that cutting plane problem, based on the sense of
the constraint, and the constant term corresponding to the deterministic part
of the constraint.
"""
function build_cut_objective_sparse(rm::Model, unc_con::UncConstraint)
    rmext      = get_robust(rm)
    unc_coeffs = Dict{Int,Float64}()
    lhs_const  = 0.0
    con_sense  = (JuMP.sense(unc_con) == :(<=)) ? :Max  :  :Min
    # Uncertain expression--variable terms
    for (unc_expr, var) in linearterms(unc_con.terms)
        var_val = getvalue(var)
        for (unc_coeff, unc) in linearterms(unc_expr)
            cur_unc_coeff = get(unc_coeffs, unc.id, 0.0)
            unc_coeffs[unc.id] = cur_unc_coeff + unc_coeff * var_val
        end
        lhs_const += unc_expr.constant * var_val
    end
    # Standalone uncertain expression/constant term
        for (unc_coeff, unc) in linearterms(unc_con.terms.constant)
            cur_unc_coeff = get(unc_coeffs, unc.id, 0.0)
            unc_coeffs[unc.id] = cur_unc_coeff + unc_coeff
        end
        lhs_const += unc_con.terms.constant.constant
    return con_sense, collect(unc_coeffs), lhs_const
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
