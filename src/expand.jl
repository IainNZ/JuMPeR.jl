#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2016: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# src/adaptive/expand.jl
# Adaptive robust optimization support - pre-solve expansion
#-----------------------------------------------------------------------

any_adaptive(u::UncVarExpr) = any(v->isa(v,Adaptive), u.vars)


function expand_adaptive(rm::Model)
    # Extract robust-specifc part of model information
    rmext = get_robust(rm)::RobustModelExt
    # If no adaptive variables, bail right away
    isempty(rmext.adp_policy) && return
    # Collect the new variables or expressions that will replace
    # the adaptive variables in the current constraints
    new_vars = Any[]
    # Collect the new constraints we'll be adding to make the
    # policies make sense
    new_cons = UncConstraint[]
    # For every adaptive variable...
    for i in 1:rmext.num_adps
        pol = rmext.adp_policy[i]  # Extract the policy
        #---------------------------------------------------------------
        if pol == :Static
            # Static policy - no dependence on the uncertain parameters
            # Create a normal variable to replace this adaptive variable
            push!(new_vars, Variable(rm,
                rmext.adp_lower[i], rmext.adp_upper[i],
                rmext.adp_cat[i], rmext.adp_names[i]) )
        #---------------------------------------------------------------
        elseif pol == :Affine
            # Extract the container of uncertain parameters that this
            # adaptive variable depends on
            deps = rmext.adp_arguments[i]
            # Ensure stage-wise dependence not being used
            #if rmext.adpStage[i] != 0
            #    error("Not implemented!")
            #end
            # Create auxiliary variables - one for each uncertain parameter,
            # and one for the independent term
            aux_aff = Dict{Any,Variable}()
            for j in eachindex(deps)
                # Construct a roughly sensible name for the auxiliary using
                # the adaptive variable and uncertain parameter names
                vname = string(adp_str(rm,i), "{", deps[j], "}")
                aux_aff[j] = Variable(rm, -Inf, +Inf, :Cont, vname)
            end
            vname = string(adp_str(rm,i), "{_}")
            aux_con = Variable(rm, -Inf, +Inf, :Cont, vname)
            # Build the policy
            aff_policy = UncExpr(1) * aux_con
            for j in eachindex(deps)
                push!(aff_policy, UncExpr(deps[j]), aux_aff[j])
            end
            #println(aff_policy)
            push!(new_vars, aff_policy)
            # Add bound constraints on the policy
            if rmext.adp_lower[i] != -Inf  # Lower bound?
                push!(new_cons, UncConstraint(aff_policy, rmext.adp_lower[i], +Inf))
            end
            if rmext.adp_upper[i] != +Inf  # Lower bound?
                push!(new_cons, UncConstraint(aff_policy, -Inf, rmext.adp_upper[i]))
            end
        #---------------------------------------------------------------
        else
            error("Unknown policy type")
        end
    end

    # Create new constraints from uncertain-and-variable constraints
    # that contained an adaptive variable
    for uncaffcon in rmext.unc_constraints
        lhs = uncaffcon.terms
        !any_adaptive(lhs) && continue
        new_lhs = UncVarExpr(lhs.constant)
        for (coeff, var) in linearterms(lhs)
            new_lhs += coeff * (isa(var, Adaptive) ? new_vars[var.id] : var)
        end
        push!(new_cons, UncConstraint(new_lhs, uncaffcon.lb, uncaffcon.ub))
        # Remove old constraint by emptying all fields
        lhs.vars = JuMPeRVar[]
        lhs.coeffs = UncExpr[]
        lhs.constant = UncExpr()
        if uncaffcon.lb != -Inf
            uncaffcon.lb = 0
        end
        if uncaffcon.ub != +Inf
            uncaffcon.ub = 0
        end
    end

    # Create new constraints from number-and-adaptive constraints
    for varaffcon in rmext.adapt_constraints
        lhs = varaffcon.terms
        new_lhs = UncVarExpr(lhs.constant)
        for (coeff, var) in linearterms(lhs)
            new_lhs += coeff * (isa(var, Adaptive) ? new_vars[var.id] : var)
        end
        push!(new_cons, UncConstraint(new_lhs, varaffcon.lb, varaffcon.ub))
        # We don't need to do anything with the old constraints
    end

    # Add these constraints to the RobustModel
    map(c->JuMP.addconstraint(rm, c), new_cons)
end
