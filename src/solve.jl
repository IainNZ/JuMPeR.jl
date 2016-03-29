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

function solve_robust(rm::Model; suppress_warnings=false, report=false,
                        active_cuts=false, add_box=false, show_cuts=false,
                        kwargs...)
    _solve_robust(rm,suppress_warnings,report,active_cuts,add_box,
                    show_cuts,kwargs)
end

function _solve_robust(rm::Model, suppress_warnings::Bool, report::Bool,
                        active_cuts::Bool, add_box::Union{Float64,Bool},
                        show_cuts::Bool, kwargs::Vector{Any})
    # Pull out extra keyword arguments that we will pas through to oracles
    prefs = Dict{Symbol,Any}([name => value for (name,value) in kwargs])
    # TEMP: Expand out adaptive terms
    expand_adaptive(rm)

    rmext = get_robust(rm)::RobustModelExt
    start_time = time()

    #-------------------------------------------------------------------
    # PART 1: MASTER SETUP
    #-------------------------------------------------------------------
    # Create master problem, a normal JuMP model
    master = Model(solver=rm.solver)
    # Variables
    master.numCols   = rm.numCols
    master.colNames  = copy(rm.colNames)
    master.colLower  = copy(rm.colLower)
    master.colUpper  = copy(rm.colUpper)
    master.colCat    = copy(rm.colCat)
    master.colVal    = copy(rm.colVal)
    mastervars       = Variable[Variable(master, i) for i in 1:rm.numCols]
    # Objective (copy, it must be certain)
    master.objSense  = rm.objSense
    master.obj       = copy_quad(rm.obj, mastervars)
    # Certain constraints
    master.linconstr = map(con -> copy_linconstr(con, mastervars), rm.linconstr)
    master.quadconstr = map(con -> copy_quadconstr(con, mastervars), rm.quadconstr)
    # Copy JuMPContainers over so we get good printing
    master.dictList  = copy(rm.dictList)
    # Copy the callbacks, even though they won't all work correctly
    if !isempty(rm.callbacks)
        if !suppress_warnings
            Base.warn("""Only informational callbacks will work as expected with JuMPeR.
                         If you wish to use lazy constraint callbacks, see the JuMPeR
                         source code for addSpecialLazyCallback() for help.""")
        end
        master.callbacks = copy(rm.callbacks)
    end
    # Create the wrapper for special lazy callbacks
    function callback_passthru(cb)
        # Copy primal solution back to outer RobustModel
        rm.colVal = copy(master.colVal)
        # Collect constraints to add, which will be for the RobustModel
        outer_cons = Any[]
        for f in rmext.lazy_callbacks
            ret = f(cb)
            ret != nothing && append!(outer_cons, ret)
        end
        # Convert constraints and add to master
        map(con ->JuMP.addLazyConstraint(cb, copy_linconstr(con, mastervars)), outer_cons)
    end
    #if !isempty(rmext.lazy_callbacks)
    #    addLazyCallback(master, callback_passthru)
    #end

    num_unccons = length(rmext.unc_constraints)

    # If the problem is a MIP, we are going to have to do more work
    is_mip = any(map(cat -> (cat in [:Int,:Bin]), master.colCat))


    # Put box around original solution.
    # Really should be doing something smarter, only need this
    # to deal with an unbounded initial solution in the cutting
    # plane process. Can cause numerical problems, so its optional.
    if add_box != false
        !(isa(add_box, Real)) && error("add_box option should be false or a number")
        map!(v -> v == -Inf ? -add_box : v, master.colLower)
        map!(v -> v ==  Inf ?  add_box : v, master.colUpper)
    end

    master_init_time = time() - start_time

    #########################################################################
    # PART 2: ORACLE SETUP

    for ind in 1:num_unccons
        # Associate constraints with the default oracle if they haven't been
        # assigned one otherwise
        if rmext.constraint_uncsets[ind] == nothing
            rmext.constraint_uncsets[ind] = rmext.default_uncset
        end
        # Register the constraints with their oracles by providing the index
        register_constraint(rmext.constraint_uncsets[ind], rm, ind, prefs)
    end

    # Give oracles time to do any setup. For example, generating the cutting
    # plane model, or determing what a reformulation should look like, or
    # processing a data set for the particular model at hand.
    oracle_setup_time = time()
    oracle_to_cons = Dict{AbstractUncertaintySet,Vector{Int}}()
    for ind in 1:num_unccons
        oracle = rmext.constraint_uncsets[ind]
        if !(oracle in keys(oracle_to_cons))
            setup_set(oracle, rm, active_cuts, prefs)
            oracle_to_cons[oracle] = [ind]
        else
            push!(oracle_to_cons[oracle], ind)
        end
    end
    oracle_setup_time = time() - oracle_setup_time

    #########################################################################
    # PART 3: REFORMULATION
    reform_time = time()
    reformed_cons = 0
    for oracle in keys(oracle_to_cons)
        reformed_cons += generate_reform(oracle, master, rm, oracle_to_cons[oracle])
    end
    reform_time = time() - reform_time

    #########################################################################
    # PART 4: CUTTING PLANES
    master_status   = nothing
    cutting_rounds  = 0
    cuts_added      = 0
    master_time     = 0.0
    cut_time        = 0.0
    #   LP: add constraints, hot starting if supported
    # MILP: lazy constraint callback
    if is_mip
        function lazyCallback(cb)
            cutting_rounds += 1
            get(prefs, :debug_printcut, false) && println("CUTTING ROUND $cutting_rounds")

            # Generate cuts
            tic()
            for oracle in keys(oracle_to_cons)
                cons_to_add = generate_cut(oracle, master, rm, oracle_to_cons[oracle])
                for new_con in cons_to_add
                    addLazyConstraint(cb, new_con)
                end
                cuts_added += length(cons_to_add)
            end
            cut_time += toq()
        end
        addLazyCallback(master, lazyCallback)

        # Solve master. Terminate when we have an optimal integer solution
        # and no lazy constraints are added
        tic()
        master_status = solve(master, suppress_warnings=true)
        master_time = toq() - cut_time
        if master_status == :Infeasible
            !suppress_warnings && Base.warn("JuMPeR: master problem (integer) is infeasible.")
        elseif master_status == :Unbounded
            !suppress_warnings && Base.warn("""JuMPeR: master problem (integer) is unbounded.
                                               This may be due to:
                                               * the problem requiring additional cutting planes for it to
                                                 become bounded. Consider adding bounds to variables
                                                 or use the add_box option to solve.
                                               * the uncertainty set being empty. Check the uncertainty
                                                 set has at least one value.
                                               * the problem actually being unbounded.""")
        end
    else # not a MIP
        # In the event of unboundedness, and if we have a ray, we will see if we get the
        # same ray twice. To do so, we will store the ray on the first iteration and
        # check against it in the next iteration.
        unbound_ray = Float64[]
        # Begin main solve loop
        while true
            cutting_rounds += 1
            get(prefs, :debug_printcut, false) && println("CUTTING ROUND $cutting_rounds")

            # Solve master
            tic()
            master_status = solve(master, suppress_warnings=true)
            master_time += toq()

            # If master is infeasible, we should stop now
            if master_status == :Infeasible
                !suppress_warnings && Base.warn("JuMPeR: master problem (continuous) is infeasible.")
                break
            end

            # If master is unbounded, we are in one of these situations:
            # 1. Problem is bounded, just needs cuts. We have a ray.
            # 2. Problem is bounded, just needs cuts. We don't have a ray.
            # 3. Problem is unbounded, we have a ray.
            # 4. Problem is unbounded, we don't have a ray.
            # First step is to look at ray situation
            if master_status == :Unbounded
                if (length(master.colVal) >0 && isnan(master.colVal[1])) ||
                   (length(master.colVal)==0)
                    # No ray available (2 or 4) - we can't add a cut even if
                    # we wanted to.
                    !suppress_warnings &&
                    Base.warn("""JuMPeR: master problem (continuous) is unbounded and
                                 no ray is available. Lack of a ray is due to the particular
                                 solver selected, but the unboundedness may be due to:
                                 * the problem requiring additional cutting planes for it to
                                   become bounded. Consider adding bounds to variables
                                   or use the add_box option to solve.
                                 * the uncertainty set being empty. Check the uncertainty
                                   set has at least one value.
                                 * the problem actually being unbounded.""")
                    break
                else
                    # We have a ray (1 or 3) - we don't know for sure what is
                    # happening. We will use a heuristic: if we get the same
                    # solution twice in a row, we will assume it is unbounded-unbounded.
                    if length(unbound_ray) == 0
                        # We haven't been unbounded yet, lets see if we get same ray
                        # twice before giving up
                        unbound_ray = copy(master.colVal)
                    else
                        # Compare with previous ray. Not clear what the tolerance
                        # should be, but lets just go with something simple until
                        # it goes haywire somewhere
                        if norm(unbound_ray .- master.colVal) <= 1e-6
                            # Same ray again
                            !suppress_warnings &&
                            Base.warn("""JuMPeR: master problem (continuous) is unbounded, but
                                         ray is available. Attempted to add cuts using this ray, but same
                                         ray was returned by solver again, so assuming problem is
                                         unbounded and terminating. Unboundedness may be due to:
                                         * the problem requiring additional cutting planes for it to
                                           become bounded. Consider adding bounds to variables
                                           or use the add_box option to solve.
                                         * the uncertainty set being empty. Check the uncertainty
                                           set has at least one value.
                                         * the problem actually being unbounded.""")
                            break
                        else
                            # Different ray, update ray in case we get same thing
                            # Note that this will fail if it cycles, but that is
                            # a really pathological case.
                            unbound_ray = copy(master.colVal)
                        end
                    end
                end
            end

            # Generate cuts
            cut_added = false
            tic()
            show_cuts && print_with_color(:white, "Cutting round $(cutting_rounds)\n")
            for oracle in keys(oracle_to_cons)
                cons_to_add = generate_cut(oracle, master, rm, oracle_to_cons[oracle])
                show_cuts && _show_cuts(rm, oracle, oracle_to_cons[oracle], cons_to_add)
                for new_con in cons_to_add
                    addConstraint(master, new_con)
                    cut_added = true
                end
                cuts_added += length(cons_to_add)
            end
            cut_time += toq()

            # Terminate solve loop when no more cuts added, with one caveat:
            # if we are unbounded, have a ray, and either
            # * using reformulation
            # * or no cuts were generated from the ray
            # then we'll never display a message from the code above.
            if master_status == :Unbounded && !cut_added
                !suppress_warnings &&
                Base.warn("""JuMPeR: master problem (continuous) is unbounded, but
                             ray is available. No cuts were added by oracles, so assuming
                             problem is unbounded and terminating. Unboundedness may be due to:
                             * the problem requiring additional cutting planes for it to
                               become bounded. Consider adding bounds to variables
                               or use the add_box option to solve.
                             * the uncertainty set being empty. Check the uncertainty
                               set has at least one value.
                             * the problem actually being unbounded.""")
            end
            !cut_added && break
        end
    end

    # Return solution
    total_time = time() - start_time
    rm.colVal = copy(master.colVal)
    rm.objVal = master.objVal

    # DEBUG: If user wants it, print final model
    if get(prefs, :debug_printfinal, false)
        println("BEGIN DEBUG :debug_printfinal")
        print(master)
        println("END DEBUG   :debug_printfinal")
    end

    # OPTION: Get active cuts (1 per constraint) ("active_cuts=true")
    # VERY EXPERIMENTAL
    tic()
    if active_cuts
        rmext.scenarios = Vector{Nullable{Scenario}}(length(rmext.unc_constraints))
        for oracle in keys(oracle_to_cons)
            scens_to_add = generate_scenario(oracle, master, rm, oracle_to_cons[oracle])
            for i in 1:length(oracle_to_cons[oracle])
                rmext.scenarios[oracle_to_cons[oracle][i]] = scens_to_add[i]
            end
        end
    end
    activecut_time = toq()

    # OPTION: Report ("report=true")
    if report
        println("Solution report")
        @printf("Uncertain cons.  %12d\n", num_unccons)
        @printf("  Reformulated   %12d\n", reformed_cons)
        @printf("Cutting rounds   %12d\n", cutting_rounds)
        @printf("Total cuts:      %12d\n", cuts_added)
        @printf("Overall time:    %12.5f\n", total_time)
        @printf("  Master init    %12.5f (%6.2f%%)\n", master_init_time, master_init_time/total_time*100)
        @printf("  Oracle setup   %12.5f (%6.2f%%)\n", oracle_setup_time, oracle_setup_time/total_time*100)
        @printf("  Reformulation  %12.5f (%6.2f%%)\n", reform_time, reform_time/total_time*100)
        @printf("  Master solve   %12.5f (%6.2f%%)\n", master_time, master_time/total_time*100)
        @printf("  Cut solve&add  %12.5f (%6.2f%%)\n", cut_time, cut_time/total_time*100)
        active_cuts && @printf("Active cut time: %12.5f\n", activecut_time)
    end

    # Store the internal model
    rm.internalModel = master.internalModel
    rm.internalModelLoaded = true

    # Return solve status
    return master_status
end


function _show_cuts(rm, oracle, con_indices, cons_added)
    rmext = get_robust(rm)
    println("Oracle: ", typeof(oracle))
    println("Constraint", length(con_indices) <= 1 ? ":" : "s:")
    for idx in con_indices
        con = rmext.unc_constraints[idx]
        println("  ", con)
    end
    if length(cons_added) == 0
        println("No generated cuts")
    else
        println("Generated cut", length(cons_added) == 1 ? ":" : "s:")
        for con in cons_added
            println("  ", con)
        end
    end
end
