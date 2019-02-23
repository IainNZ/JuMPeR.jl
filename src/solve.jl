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

function solve_robust(rm::Model; suppress_warnings=false,
                        disable_cuts=false, active_scenarios=false,
                        show_cuts=false, kwargs...)
    _solve_robust(rm,suppress_warnings,disable_cuts,active_scenarios,show_cuts,kwargs)
end

function _solve_robust(rm::Model, suppress_warnings::Bool,
                        disable_cuts::Bool, active_scenarios::Bool,
                        show_cuts::Bool, kwargs)

    rmext = get_robust(rm)::RobustModelExt
    if rmext.solved
        error("JuMPeR: RobustModel has already been solved. Unlike JuMP Models,
               RobustModels are permanently modified when solve is called. If
               you need to solve a succession of robust optimization problems,
               you should build the model from scratch each time.")
    end

    #-------------------------------------------------------------------
    # 1. Expand out adaptive terms
    #-------------------------------------------------------------------
    expand_adaptive(rm)

    #-------------------------------------------------------------------
    # 2. Setup for uncertainty sets
    #-------------------------------------------------------------------
    # Prepare to pass through preferences to a sets
    prefs = Dict{Symbol,Any}(name => value for (name,value) in kwargs)
    # Build mapping from uncertainty sets to constraints
    uncsets_to_con_idxs = Dict{AbstractUncertaintySet,Vector{Int}}()
    for idx in 1:length(rmext.constraint_uncsets)
        # Associate constraints with the default set if they haven't
        # been assigned one manually
        cur_uncset = rmext.constraint_uncsets[idx]
        if cur_uncset == nothing
            rmext.constraint_uncsets[idx] = rmext.default_uncset
            cur_uncset = rmext.default_uncset
        end
        if cur_uncset in keys(uncsets_to_con_idxs)
            push!(uncsets_to_con_idxs[cur_uncset], idx)
        else
            uncsets_to_con_idxs[cur_uncset] = [idx]
        end
    end
    # Give uncertainty sets chances to do any pre-solve setup
    for (uncset, idxs) in uncsets_to_con_idxs
        setup_set(uncset, rm, idxs, active_scenarios, prefs)
    end

    #-------------------------------------------------------------------
    # 3. Reformulations
    #-------------------------------------------------------------------
    for (uncset, idxs) in uncsets_to_con_idxs
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
                for (uncset, idxs) in uncsets_to_con_idxs
                    cons_to_add = generate_cut(uncset, rm, idxs)
                    for new_con in cons_to_add
                        JuMP.addlazyconstraint(cb, new_con)
                    end
                end
            end
            addlazycallback(rm, lazyCallback)
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
            for (uncset, idxs) in uncsets_to_con_idxs
                cons_to_add = generate_cut(uncset, rm, idxs)
                for new_con in cons_to_add
                    JuMP.addconstraint(rm, new_con)
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
    # 5. Extract active scenarios, if desired
    #-------------------------------------------------------------------
    if active_scenarios
        rmext.scenarios = Vector{Union{Scenario, Missing}}(undef, length(rmext.unc_constraints))
        for (uncset, idxs) in uncsets_to_con_idxs
            scens_to_add = generate_scenario(uncset, rm, idxs)
            for i in 1:length(idxs)
                rmext.scenarios[idxs[i]] = scens_to_add[i]
            end
        end
    end

    # Return solve status
    return status
end
