#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2016: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# src/solve.jl
# Core of the logic for solving robust optimization problems. Uses
# uncertainty sets to build up the deterministic equivalent problem.
# Included by src/JuMPeR.jl
#-----------------------------------------------------------------------

# Warning messages
include("solve_msgs.jl")

# Utility functions for moving an expression to a new model
copy_quad(old_quad::QuadExpr, new_vars) =
    QuadExpr(   Variable[new_vars[v.col] for v in old_quad.qvars1],
                Variable[new_vars[v.col] for v in old_quad.qvars2],
                copy(old_quad.qcoeffs),
                copy_aff(old_quad.aff, new_vars))
copy_aff(old_aff::AffExpr, new_vars) =
    AffExpr(    Variable[new_vars[v.col] for v in old_aff.vars],
                copy(old_aff.coeffs),
                old_aff.constant)
copy_linconstr(old_con::LinearConstraint, new_vars) =
    LinearConstraint(copy_aff(old_con.terms, new_vars), old_con.lb, old_con.ub)
copy_quadconstr(old_con::QuadConstraint, new_vars) =
    QuadConstraint(copy_quad(old_con.terms, new_vars), old_con.sense)

#-----------------------------------------------------------------------

function solve_robust(rm::Model; suppress_warnings=false,
                        disable_cuts=false, active_scenarios=false,
                        show_cuts=false, kwargs...)
    _solve_robust(rm,suppress_warnings,disable_cuts,active_scenarios,show_cuts,kwargs)
end

function _solve_robust(rm::Model, suppress_warnings::Bool,
                        disable_cuts::Bool, active_scenarios::Bool,
                        show_cuts::Bool, kwargs::Vector{Any})

    rmext = get_robust(rm)::RobustModelExt
    if rmext.solved
        error("JuMPeR: RobustModel has already been solved. Unlike JuMP Models,
               RobustModels are permanently modified when solve is called. If
               you need to solve a succession of robust optimization problems,
               you should build the model from scratch each time.")
    end

    # TEMP: Expand out adaptive terms
    expand_adaptive(rm)

    #-------------------------------------------------------------------
    # 2. Setup for uncertainty sets
    #-------------------------------------------------------------------
    # Prepare to pass through preferences to a sets
    prefs = Dict{Symbol,Any}([name => value for (name,value) in kwargs])
    # Register the constraints with their uncertainty sets
    for idx in 1:length(rmext.constraint_uncsets)
        # Associate constraints with the default set if they haven't
        # been assigned one manually
        if rmext.constraint_uncsets[idx] == nothing
            rmext.constraint_uncsets[idx] = rmext.default_uncset
        end
        register_constraint(rmext.constraint_uncsets[idx], rm, idx, prefs)
    end
    # Give uncertainty sets chances to do any pre-solve setup
    uncsets_to_cons = Dict{AbstractUncertaintySet,Vector{Int}}()
    for idx in 1:length(rmext.constraint_uncsets)
        uncset = rmext.constraint_uncsets[idx]
        if uncset ∉ keys(uncsets_to_cons)
            # Uncset hasn't been set up yet
            setup_set(uncset, rm, active_scenarios, prefs)
            uncsets_to_cons[uncset] = [idx]
        else
            # Already done set up
            push!(uncsets_to_cons[uncset], idx)
        end
    end

    #-------------------------------------------------------------------
    # 3. Reformulations
    #-------------------------------------------------------------------
    for (uncset, idxs) in uncsets_to_cons
        generate_reform(uncset, rm, idxs)
    end

    #-------------------------------------------------------------------
    # 4. Cutting planes and main solve loop
    #   LP: add constraints, hot starting if supported
    # MILP: lazy constraint callback
    #-------------------------------------------------------------------
    status = :Unsolved
    cutting_rounds = 0
    # If the problem is a MIO...
    if any(map(cat -> (cat in [:Int,:Bin]), rm.colCat))
        # Unless cutting planes have been explictly disabled...
        if !disable_cuts
            # ... add a cutting plane callback
            function lazyCallback(cb)
                cutting_rounds += 1
                show_cuts && println("JuMPeR: Cutting plane callback $cutting_rounds")
                for (uncset, idxs) in uncsets_to_cons
                    cons_to_add = generate_cut(uncset, rm, idxs)
                    for new_con in cons_to_add
                        addLazyConstraint(cb, new_con)
                    end
                end
            end
            addLazyCallback(rm, lazyCallback)
        end
        # Solve the problem, terminating when we have an optimal integer
        # solution and no lazy constraints are added. Be prepared for the case
        # where the solver throws an error, probably due to not supporting
        # lazy contraints (e.g. Cbc)
        status = :Error
        try
            status = solve(rm, ignore_solve_hook=true, suppress_warnings=true)
        catch err
            if !suppress_warnings
                warn_ip_error()
                println(err)
            end
        end
        if !suppress_warnings
            status == :Infeasible && warn_ip_infeasible()
            status == :Unbounded  && warn_ip_unbounded()
        end
    # If the problem is not a MIO...
    else
        # In the event of unboundedness, and if we have a ray, we will see if
        # we get the same ray twice. To do so, we will store the ray on the
        # first iteration and check against it in the next iteration.
        unbound_ray = Float64[]
        # Begin main solve loop
        while true
            cutting_rounds += 1
            show_cuts && println("JuMPeR: Cutting plane round $cutting_rounds")
            # Solve the current relaxation
            status = solve(rm, ignore_solve_hook=true, suppress_warnings=true)
            # If it is already infeasible, we should stop now
            if status == :Infeasible
                !suppress_warnings && warn_lp_infeasible()
                break
            end
            # If it is unbounded, we are in one of these situations:
            # 1. Problem is bounded, just needs cuts. We have a ray.
            # 2. Problem is bounded, just needs cuts. We don't have a ray.
            # 3. Problem is unbounded, we have a ray.
            # 4. Problem is unbounded, we don't have a ray.
            # 5. Problem is unbounded and we disabled cuts, so stop.
            # First, determine ray situation
            if status == :Unbounded
                # If we're in case 5, we can just stop
                if disable_cuts
                    !suppress_warnings && warn_lp_unbounded()
                    break
                end
                # If we have decision variable values but they are NaN, or if
                # we have no decision variable values...
                if (length(rm.colVal) >0 && isnan(rm.colVal[1])) ||
                   (length(rm.colVal)==0)
                    # No ray available (cases 2 or 4), so we couldn't add a cut
                    # even if we wanted to due to lack of information.
                    !suppress_warnings && warn_lp_unbounded_noray()
                    break
                else
                    # We have a ray (cases 1 or 3), so we don't know for sure
                    # what is happening (really is an unbounded problem, or
                    # just need some more cutting planes). We will use a
                    # heuristic: if we get the same solution twice in a row,
                    # we will assume it is unbounded-unbounded.
                    if length(unbound_ray) == 0
                        # We haven't been unbounded yet, lets see if we get
                        # same ray twice before giving up
                        unbound_ray = copy(rm.colVal)
                    else
                        # Compare with previous ray. Not clear what the tolerance
                        # should be, but just go with something simple until
                        # there is a good reason to change it
                        if norm(unbound_ray - rm.colVal) <= 1e-6
                            # Same ray again
                            !suppress_warnings && warn_lp_unbounded_sameray()
                            break
                        else
                            # Different ray, update ray in case we get same thing
                            # Note that this will fail if it cycles, but that is
                            # a really pathological case.
                            unbound_ray = copy(rm.colVal)
                        end
                    end
                end
            end  # if status == :Unbounded
            # If no cuts wanted, we can stop right here
            disable_cuts && break
            # Attempt to generate cut
            cut_added = false
            for (uncset, idxs) in uncsets_to_cons
                cons_to_add = generate_cut(uncset, rm, idxs)
                for new_con in cons_to_add
                    addConstraint(rm, new_con)
                    cut_added = true
                end
            end
            # Terminate solve loop when no more cuts added, with one caveat:
            # if we are unbounded, have a ray, and either
            # * using reformulation
            # * or no cuts were generated from the ray
            # then we'll never display any message from the code above.
            if status == :Unbounded && !cut_added && !suppress_warnings
                warn_lp_unbounded()
            end
            # Stop once no more cuts are added by the sets
            !cut_added && break
        end
    end

    #-------------------------------------------------------------------
    # 4. Extract active scenarios, if desired
    #-------------------------------------------------------------------
    if active_scenarios
        rmext.scenarios = Vector{Nullable{Scenario}}(length(rmext.unc_constraints))
        for (uncset, idxs) in uncsets_to_cons
            scens_to_add = generate_scenario(uncset, rm, idxs)
            for i in 1:length(idxs)
                rmext.scenarios[idxs[i]] = scens_to_add[i]
            end
        end
    end

    # Return solve status
    return status
end
