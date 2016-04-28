#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2016: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# src/uncsets_basic_cuts.jl
# Cutting plane-specific functionality for the BasicUncertaintySet,
# including setup_set_cut and generate_cut.
# Included by src/uncsets_basic.jl
#-----------------------------------------------------------------------


"""
    setup_set_cut(BasicUncertaintySet, RobustModel)

Builds the cutting plane model for the BasicUncertaintySet.
"""
function setup_set_cut(us::BasicUncertaintySet, rm::Model)
    # The uncertain parameters of the RobustModel become the variables
    # of the cutting plane model.
    rmext = get_robust(rm)
    cut_model = Model()
    if isa(rmext.cutsolver, JuMP.UnsetSolver)
        cut_model.solver = rm.solver
    else
        cut_model.solver = rmext.cutsolver
    end
    cut_model.numCols  = copy(rmext.num_uncs)
    cut_model.colNames = copy(rmext.unc_names)
    cut_model.colLower = copy(rmext.unc_lower)
    cut_model.colUpper = copy(rmext.unc_upper)
    cut_model.colCat   = copy(rmext.unc_cat)

    # Create a vector of Variable objects for easier translation of
    # the uncertainty set constraints into deterministic constraints
    cut_vars = [Variable(cut_model, i) for i in 1:rmext.num_uncs]

    # Convert linear constraints
    for c in us.linear_constraints
        lhs = uaff_to_aff(c.terms, cut_vars)
        new_c = LinearConstraint(lhs, c.lb, c.ub)
        JuMP.addconstraint(cut_model, new_c)
    end

    # Convert norm constraints
    for (idx, norm_c) in enumerate(us.norm_constraints)
        # Norm constraints
        # Input:  ‖[a₁ᵀu, a₂ᵀu, ...]‖ ≤ Γ
        normexp = norm_c.normexpr
        terms   = normexp.norm.terms  # Inside of ‖…‖
        Γ       = -normexp.aff.constant / normexp.coeff
        n_terms = length(terms)
        if isa(norm_c, UncSetNormConstraint{2})
            # L2 norm constraint
            # Output: yᵢ = aᵢᵀu, t = Γ, Σyᵢ^2 ≤ t^2
            @variable(cut_model, y[i=1:n_terms], basename="_c$(idx)_L2_y")
            for i in 1:n_terms
                @constraint(cut_model, y[i] == uaff_to_aff(terms[i], cut_vars))
            end
            @variable(cut_model, t == Γ, basename="_c$(idx)_L2_t")
            @constraint(cut_model, dot(y,y) <= t^2)
        elseif isa(norm_c, UncSetNormConstraint{1})
            # L1 norm constraint
            # Output: yᵢ ≥ aᵢᵀu, yᵢ ≥ -aᵢᵀu, ∑yᵢ ≤ Γ
            @variable(cut_model, y[i=1:n_terms], basename="_c$(idx)_L1_y")
            for i in 1:n_terms
                rhs = uaff_to_aff(terms[i], cut_vars)
                @constraint(cut_model, y[i] ≥ +rhs)
                @constraint(cut_model, y[i] ≥ -rhs)
            end
            @constraint(cut_model, sum(y) <= Γ)
        elseif isa(norm_c, UncSetNormConstraint{Inf})
            # L∞ norm constraint
            # Output: aᵢᵀu ≤ Γ, aᵢᵀu ≥ -Γ
            for i in 1:n_terms
                lhs = uaff_to_aff(terms[i], cut_vars)
                @constraint(cut_model, lhs ≤ +Γ)
                @constraint(cut_model, lhs ≥ -Γ)
            end
        else
            error("Unrecognized norm in uncertainty set!")
        end
    end  # enumerate(us.norm_constraints)

    us.cut_model = cut_model
    us.cut_vars = cut_vars
end


"""
    get_worst_case_value(BasicUncertaintySet, ...)

Internal function. For a given constraint and solution, find the worst-case
values of the uncertain parameters. This function is called by both
`generate_scenario` and `generate_cut`
"""
function get_worst_case_value(us::BasicUncertaintySet, rm::Model, idx::Int)
    rmext = get_robust(rm)::RobustModelExt
    con = rmext.unc_constraints[idx]
    # Update the cutting plane problem's objective, and solve
    cut_sense, unc_obj_coeffs, lhs_const = JuMPeR.build_cut_objective_sparse(rm, con)
    @objective(us.cut_model, cut_sense, sum{u[2]*us.cut_vars[u[1]], u=unc_obj_coeffs})
    cut_solve_status = solve(us.cut_model, suppress_warnings=true)
    cut_solve_status != :Optimal &&
        error("BasicUncertaintySet: cutting plane problem is infeasible or unbounded!")
    lhs_of_cut = getobjectivevalue(us.cut_model) + lhs_const

    return lhs_of_cut, us.cut_model.colVal
end


"""
    generate_scenario(BasicUncertaintySet, ...)

Wraps the results from `get_worst_case_value` in `Scenario` objects.
"""
function generate_scenario(us::BasicUncertaintySet, rm::Model, idxs::Vector{Int})
    # We need to return one Scenario per constraint
    scens = Nullable{Scenario}[]
    for idx in idxs
        _, uncvalues = get_worst_case_value(us, rm, idx)
        scen = Scenario(uncvalues)
        push!(scens, Nullable{Scenario}(scen))
    end
    return scens
end

"""
    generate_cut(BasicUncertaintySet, ...)

Given constraints with uncertain parameters, create new constraints if they
would cause the current solution to be infeasible.
"""
function generate_cut(us::BasicUncertaintySet, rm::Model, idxs::Vector{Int})
    # If not doing cuts, then fail fast
    if !us.use_cuts
        return Any[]
    end
    rmext = get_robust(rm)
    new_cons = Any[]
    for idx in idxs
        con = rmext.unc_constraints[idx]
        lhs_of_cut, uncvalues = get_worst_case_value(us, rm, idx)
        # Check violation
        if check_cut_status(con, lhs_of_cut, us.cut_tol) != :Violate
            continue  # No violation, no new cut
        end
        # Create and add the new constraint
        new_con = JuMPeR.build_certain_constraint(con, uncvalues)
        push!(new_cons, new_con)
    end
    return new_cons
end
