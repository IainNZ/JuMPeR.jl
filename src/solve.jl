#############################################################################
# JuMPeR
# Julia for Mathematical Programming - extension for Robust Optimization
# See http://github.com/IainNZ/JuMPeR.jl
#############################################################################
# solve.jl
# All the logic for solving robust optimization problems. Communicates with
# oracles to build and solve the robust counterpart.
#############################################################################


convert_model!(old_con::GenericRangeConstraint, new_m::Model) =
    map((v)->(v.m = new_m), old_con.terms.vars)
function convert_model!(old_obj::QuadExpr, new_m::Model)
    map((v)->(v.m = new_m), old_obj.qvars1)
    map((v)->(v.m = new_m), old_obj.qvars2)
    map((v)->(v.m = new_m), old_obj.aff.vars)
end


function solveRobust(rm::Model; report=false, active_cuts=false, args...)
    
    prefs = Dict()
    for (name,value) in args
        prefs[name] = value
    end

    if get(prefs, :debug_printfinal, false) ||
       get(prefs, :debug_printcut, false)
       JuMP.fillVarNames(rm)
       JuMPeR.fillUncNames(rm)
    end

    robdata = getRobust(rm)
    
    ###########################################################################
    # Expand adaptability in place
    # First, create the auxiliary expressions
    aux_expr = Dict{Int,FullAffExpr}()
    for col in keys(robdata.adapt_type)
        new_faff = FullAffExpr()
        const_term = Variable(rm, -Inf, +Inf, 0, rm.colNames[col]*"_const")
        new_faff += const_term
        for unc in robdata.adapt_on[col]
            unc_term = Variable(rm, -Inf, +Inf, 0, rm.colNames[col]*"_"*robdata.uncNames[unc.unc])
            new_faff += unc * unc_term
        end
        aux_expr[col] = new_faff
        rm.colLower[col] != -Inf  &&  addConstraint(rm, aux_expr[col] >= rm.colLower[col])
        rm.colUpper[col] != +Inf  &&  addConstraint(rm, aux_expr[col] <= rm.colUpper[col])
    end
    # Make replacements in constraints that already have uncertainties
    for unc_con in robdata.uncertainconstr
        new_terms = FullAffExpr()
        for ind in 1:length(unc_con.terms.vars)
            col = unc_con.terms.vars[ind].col
            if col in keys(aux_expr)
                # Replace it
                new_terms += unc_con.terms.coeffs[ind].constant * aux_expr[col]
            else
                # Leave it
                new_terms += unc_con.terms.coeffs[ind] * unc_con.terms.vars[ind]
            end
        end
        new_terms += unc_con.terms.constant
        unc_con.terms = new_terms
    end
    # Make replacements in "certain" constraints
    remove_cons = Int[]
    for con_ind in 1:length(rm.linconstr)
        cert_con = rm.linconstr[con_ind]
        # Just check first - don't build replacement unless we need it
        needs_replacements = false
        for ind in 1:length(cert_con.terms.vars)
            col = cert_con.terms.vars[ind].col
            if col in keys(aux_expr)
                needs_replacements = true
                break
            end
        end
        !needs_replacements && continue
        # We need it
        new_terms = FullAffExpr()
        for ind in 1:length(cert_con.terms.vars)
            col = cert_con.terms.vars[ind].col
            if col in keys(aux_expr)
                new_terms += cert_con.terms.coeffs[ind] * aux_expr[col]                
            else
                new_terms += cert_con.terms.coeffs[ind] * cert_con.terms.vars[ind]
            end
        end
        new_terms += cert_con.terms.constant
        push!(rm.ext[:Robust].uncertainconstr, UncConstraint(new_terms, cert_con.lb, cert_con.ub))
        push!(rm.ext[:Robust].oracles, nothing)
        push!(remove_cons, con_ind)
    end
    # Remove certain constraints which have been robustified
    # TODO: deleteat! in Julia v0.3?
    robdata = getRobust(rm)
    for ind = length(remove_cons):-1:1
        splice!(rm.linconstr, remove_cons[ind])
    end

    # Create master problem, a normal JuMP model
    # TODO: make a copy of certain things
    # TODO: Do we even need to make a copy anymore?
    start_time = time()
    master = Model(solver=rm.solver)
    master.objSense  = rm.objSense
    master.obj       = rm.obj  # For now, only certain obj
    master.linconstr = rm.linconstr
    master.numCols   = rm.numCols
    master.colNames  = rm.colNames
    master.colLower  = rm.colLower
    master.colUpper  = rm.colUpper
    master.colCat    = rm.colCat
    master.colVal    = copy(rm.colVal)
    mastervars       = [Variable(master, i) for i = 1:rm.numCols]
    master_init_time = time() - start_time
    num_unccons      = length(robdata.uncertainconstr)
    convert_model!(master.obj, master)
    for c in master.linconstr
        convert_model!(c, master)
    end

    # If the problem is a MIP, we are going to have to do more work
    isIP = false
    for cat in master.colCat
        if cat == JuMP.INTEGER
            isIP = true
            break
        end
    end

    # Add constraints based on the provided scenarios
    for scen in robdata.scenarios
        for ind in 1:num_unccons
            add_scenario_master(master, 
                                robdata.uncertainconstr[ind],
                                scen)
        end
    end

    # As a more general question, need to figure out a principled way of putting
    # box around original solution, or doing something when original solution is unbounded.
    for j in 1:master.numCols
        master.colLower[j] = max(master.colLower[j], -1e6)
        master.colUpper[j] = min(master.colUpper[j], +1e6)
    end

    # Some constraints may not have oracles, so we will create a default
    # PolyhedralOracle that is shared amongst them.
    # Register the constraints with their oracles
    oracle_register_time = time()
    for ind in 1:num_unccons
        c = robdata.uncertainconstr[ind]
        w = robdata.oracles[ind]
        if w == nothing
            w = robdata.defaultOracle
            robdata.oracles[ind] = w
        end
        registerConstraint(w, c, ind, prefs)
    end
    oracle_register_time = time() - oracle_register_time


    # Give oracles time to do any setup
    oracle_setup_time = time()
    for ind in 1:num_unccons
        setup(robdata.oracles[ind], rm)
    end
    oracle_setup_time = time() - oracle_setup_time


    # For oracles that want/have to reformulate, process them now
    reform_time = time()
    reformed_cons = 0
    for ind = 1:num_unccons
        if generateReform(robdata.oracles[ind], rm, ind, master)
            reformed_cons += 1
        end
    end
    reform_time = time() - reform_time

    #########################################################################
    # Main solve loop
    master_status   = nothing
    cutting_rounds  = 0
    cuts_added      = 0
    master_time     = 0
    cut_time        = 0
    ever_cut        = zeros(Bool, num_unccons)
    # No integers: add cuts, use warm starting
    # Intgers: user lazy constraints
    if isIP
        function lazyCallback(cb)
            cutting_rounds += 1

            # Master is solved

            # Generate cuts
            tic()
            for ind = 1:num_unccons
                num_cuts_added = generateCut(robdata.oracles[ind], rm, ind, master, cb)
                if num_cuts_added > 0
                    ever_cut[ind] = true
                    cuts_added += num_cuts_added
                end
            end
            cut_time += toq()

            # Solve will automatically terminate when we finish solve
            # and no lazy constraint are added
        end
        setLazyCallback(master, lazyCallback)

        # Solve master (timing will be for whole solve time, but we'll subtract
        # cut time to approximate)
        tic()
        master_status = solve(master,
                              load_model_only   = get(prefs, :load_model_only, false),
                              suppress_warnings = get(prefs, :suppress_warnings, false))
        master_time += toq()

        master_time -= cut_time

    else
        # Begin main solve loop
        while true
            cutting_rounds += 1
            if get(prefs, :debug_printcut, false)
                println("CUTTING ROUND $cutting_rounds")
            end
            
            # Solve master
            tic()
            master_status = solve(master,
                                  load_model_only   = get(prefs, :load_model_only, false),
                                  suppress_warnings = get(prefs, :suppress_warnings, false))
            master_time += toq()

            # Generate cuts
            cut_added = false
            tic()
            for ind = 1:num_unccons
                num_cuts_added = generateCut(robdata.oracles[ind], rm, ind, master)
                if num_cuts_added > 0
                    cut_added = true
                    cuts_added += num_cuts_added
                end
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
    tic()
    if active_cuts
        for ind = 1:num_unccons
            generateCut(robdata.oracles[ind], rm, ind, master, nothing, true)
        end
    end
    activecut_time = toq()

    # OPTION: Report ("report=true")
    if report
        println("Solution report")
        #println("Prefered method: $(preferred_mode==:Cut ? "Cuts" : "Reformulations")")
        println("Uncertain Constraints:")
        println("  No. constraints  $num_unccons")
        println("  Reformulated     $reformed_cons")
        count_cut = sum([ever_cut[i] ? 1 : 0 for i in 1:num_unccons])
        println("  Cutting plane    $count_cut")
        if isIP
            println("Lazy callbacks:  $cutting_rounds")
        else
            println("Cutting rounds:  $cutting_rounds")
        end
        println("Total cuts:      $cuts_added")
        println("Overall time:    $total_time")
        println("  Master init      $master_init_time")
        println("  Oracle setup     $oracle_setup_time")
        println("  Reformulation    $reform_time")
        if isIP
            println("  Master solve     $master_time (** for IP, == TotalMIPTime - CutTime)")
        else
            println("  Master solve     $master_time")
        end
        println("  Cut solve&add    $cut_time")
        active_cuts && println("Active cuts time:  $activecut_time")
    end

    # Return solve status
    return master_status
end


#############################################################################

function add_scenario_master(master::Model, con::UncConstraint, scenario::Dict)
    new_lhs = AffExpr()
    for term_index in 1:length(con.terms.vars)
        new_coeff = con.terms.coeffs[term_index].constant
        for unc_index in 1:length(con.terms.coeffs[term_index].vars)
            unc = con.terms.coeffs[term_index].vars[unc_index]
            if unc in keys(scenario)
                new_coeff += scenario[unc]*con.terms.coeffs[term_index].coeffs[unc_index]
            else
                println("term ", unc)
                return false
            end
        end
        push!(new_lhs, new_coeff, con.terms.vars[term_index])
    end
        new_constant = 0.0
        for unc_index in 1:length(con.terms.constant.vars)
            unc = con.terms.constant.vars[unc_index]
            if unc in keys(scenario)
                new_constant += scenario[unc]*con.terms.constant.coeffs[unc_index]
            else
                return false
            end
        end
    if sense(con) == :<=
        addConstraint(master, new_lhs + new_constant <= con.ub)
    else
        addConstraint(master, new_lhs + new_constant >= con.lb)
    end
    return true
end
