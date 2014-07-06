#############################################################################
# JuMPeR
# Julia for Mathematical Programming - extension for Robust Optimization
# See http://github.com/IainNZ/JuMPeR.jl
#############################################################################
# solve.jl
# All the logic for solving robust optimization problems. Communicates with
# oracles to build and solve the robust counterpart.
#############################################################################

const BOXSIZE = 1e6

# Utility functions for 
convert_model!(old_con::GenericRangeConstraint, new_m::Model) =
    map((v)->(v.m = new_m), old_con.terms.vars)
function convert_model!(old_obj::QuadExpr, new_m::Model)
    map((v)->(v.m = new_m), old_obj.qvars1)
    map((v)->(v.m = new_m), old_obj.qvars2)
    map((v)->(v.m = new_m), old_obj.aff.vars)
end

copy_quad(old_quad::QuadExpr, new_vars) =
    QuadExpr(   map(v -> new_vars[v.col], old_quad.qvars1),
                map(v -> new_vars[v.col], old_quad.qvars2),
                copy(old_quad.qcoeffs),
                copy_aff(old_quad.aff, new_vars))
copy_aff(old_aff::AffExpr, new_vars) =
    AffExpr(    map(v -> new_vars[v.col], old_aff.vars),
                copy(old_aff.coeffs),
                old_aff.constant)
copy_linconstr(old_con::LinearConstraint, new_vars) =
    LinearConstraint(copy_aff(old_con.terms, new_vars), old_con.lb, old_con.ub)

#############################################################################

function solveRobust(rm::Model; report=false, active_cuts=false, kwargs...)
    robdata = getRobust(rm)

    # Pull out extra keyword arguments that we will pas through to oracles
    prefs = [name => value for (name,value) in kwargs]

    # If we are doing any printing, we will need variables names. This ensures
    # that the variable name vectors get filled up
    if get(prefs, :debug_printfinal, false) || get(prefs, :debug_printcut, false)
       JuMP.fillVarNames(rm)
       JuMPeR.fillUncNames(rm)
    end

    start_time = time()

    #########################################################################
    # PART 1: MASTER SETUP
    # Create master problem, a normal JuMP model
    master = Model(solver=rm.solver)
    # Variables
    master.numCols   = rm.numCols
    master.colNames  = copy(rm.colNames)
    master.colLower  = copy(rm.colLower)
    master.colUpper  = copy(rm.colUpper)
    master.colCat    = copy(rm.colCat)
    master.colVal    = copy(rm.colVal)
    mastervars       = Variable[Variable(master, i) for i = 1:rm.numCols]
    # Objective (copy, it must be certain)
    master.objSense  = rm.objSense
    master.obj       = copy_quad(rm.obj, mastervars)
    # Certain constraints
    master.linconstr = map(con -> copy_linconstr(con, mastervars), rm.linconstr)
   
    num_unccons      = length(robdata.uncertainconstr)

    # If the problem is a MIP, we are going to have to do more work
    is_mip = any(map(cat -> cat == JuMP.INTEGER, master.colCat))

    # Add constraints based on the provided scenarios
    for scen in robdata.scenarios
        for ind in 1:num_unccons
            con = robdata.uncertainconstr[ind]
            !scen_satisfies_con(scen, con) && continue
            addConstraint(master, 
                copy_linconstr(
                    build_certain_constraint(con, scen_to_vec(scen)),
                mastervars))
        end
    end

    # Put box around original solution. Really should be doing something
    # smarter, only need this to deal with an unbounded initial solution.
    # Alternatively, relax this later in solve process
    map!(v -> v == -Inf ? -BOXSIZE : v, master.colLower)
    map!(v -> v ==  Inf ?  BOXSIZE : v, master.colUpper)

    master_init_time = time() - start_time

    #########################################################################
    # PART 2: ORACLE SETUP
    
    for ind in 1:num_unccons
        # Associate constraints with the default oracle if they haven't been
        # assigned one otherwise
        if robdata.oracles[ind] == nothing
            robdata.oracles[ind] = robdata.defaultOracle
        end
        # Register the constraints with their oracles by providing the index
        registerConstraint(robdata.oracles[ind], rm, ind, prefs)
    end

    # Give oracles time to do any setup. For example, generating the cutting
    # plane model, or determing what a reformulation should look like, or
    # processing a data set for the particular model at hand.
    oracle_setup_time = time()
    oracle_to_cons = Dict{AbstractOracle,Vector{Int}}()
    for ind in 1:num_unccons
        oracle = robdata.oracles[ind]
        if !(oracle in keys(oracle_to_cons))
            setup(oracle, rm, prefs)
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
        reformed_cons += generateReform(oracle, master, rm, oracle_to_cons[oracle])
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
                cons_to_add = generateCut(oracle, master, rm, oracle_to_cons[oracle])
                for new_con in cons_to_add
                    convert_model!(new_con, master)
                    addLazyConstraint(cb, new_con)
                end
                cuts_added += length(cons_to_add)
            end
            cut_time += toq()
        end
        setLazyCallback(master, lazyCallback)

        # Solve master. Terminate when we have an optimal integer solution
        # and no lazy constraints are added
        tic()
        master_status = solve(master,
                              load_model_only   = get(prefs, :load_model_only, false),
                              suppress_warnings = get(prefs, :suppress_warnings, false))
        master_time = toq() - cut_time
    else
        # Begin main solve loop
        while true
            cutting_rounds += 1
            get(prefs, :debug_printcut, false) && println("CUTTING ROUND $cutting_rounds")
            
            # Solve master
            tic()
            master_status = solve(master,
                                  load_model_only   = get(prefs, :load_model_only, false),
                                  suppress_warnings = get(prefs, :suppress_warnings, false))
            master_time += toq()

            # Generate cuts
            cut_added = false
            tic()
            for oracle in keys(oracle_to_cons)
                cons_to_add = generateCut(oracle, master, rm, oracle_to_cons[oracle])
                for new_con in cons_to_add
                    convert_model!(new_con, master)
                    addConstraint(master, new_con)
                    cut_added = true
                end
                cuts_added += length(cons_to_add)
            end
            cut_time += toq()

            # Terminate solve loop when no more cuts added
            if !cut_added
                break
            end
        end
    end

    # Return solution
    total_time = time() - start_time
    rm.colVal = master.colVal[1:rm.numCols]
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
        for oracle in keys(oracle_to_cons)
            generateCut(oracle, master, rm, oracle_to_cons[oracle], true)
        end
    end
    activecut_time = toq()

    # OPTION: Report ("report=true")
    if report
        println("Solution report")
        println("Uncertain cons.    $num_unccons")
        println("  Reformulated     $reformed_cons")
        println("Cutting rounds     $cutting_rounds")
        println("Total cuts:        $cuts_added")
        println("Overall time:      $total_time")
        println("  Master init      $master_init_time")
        println("  Oracle setup     $oracle_setup_time")
        println("  Reformulation    $reform_time")
        println("  Master solve     $master_time")
        println("  Cut solve&add    $cut_time")
        active_cuts && println("Active cuts time:  $activecut_time")
    end

    # Return solve status
    return master_status
end