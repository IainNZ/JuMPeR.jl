#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2015: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# src/oracle_gen.jl
# The "default" oracle - using the bounds and constraints on the
# uncertain parameters, either create and repeated solve the cutting
# plane model, or refomulate all the constraints using duality.
#-----------------------------------------------------------------------

type GeneralOracle <: AbstractOracle
    use_cuts::Bool

    # Cutting plane-specific options
    cut_model::Model
    cut_vars::Vector{Variable}
    cut_tol::Float64

    # Reformulation structure, generated in setup
    num_dualvar::Int
    dual_A::Vector{Vector{@compat Tuple{Int,Float64}}}
    dual_objs::Vector{Float64}
    dual_vartype::Vector{Symbol}
    dual_contype::Vector{Symbol}
    dual_ell_rhs_idxs::Vector{Vector{Int}}
    dual_ell_lhs_idxs::Vector{Int}
    dual_l1_rhs_idxs::Vector{Vector{Int}}
    dual_l1_lhs_idxs::Vector{Int}
    dual_ω_idxs::Vector{Vector{Int}}
    dual_ω′_idxs::Vector{Vector{Int}}

    # Options
    debug_printcut::Bool
end
# Default constructor
GeneralOracle() = 
    GeneralOracle(  false,
                    Model(), Variable[], 0.0, # Cutting plane
                    0, Vector{@compat Tuple{Int,Float64}}[], Float64[], Symbol[], Symbol[], 
                    Vector{Int}[], Int[],           # Reformulation, 2-norm
                    Vector{Int}[], Int[],           # Reformulation, 1-norm
                    Vector{Int}[], Vector{Int}[],   # Reformulation, ∞-norm
                    false)


# registerConstraint
# We must handle this constraint, and the users preferences have been
# communicated through prefs. We don't need to take any action here.
function registerConstraint(gen::GeneralOracle, rm::Model,
                                ind::Int, prefs::Dict{Symbol,Any})
    nothing
end

#-----------------------------------------------------------------------
# setup
# Generate the cutting plane model or precompute the reformulation's structure.
function setup(gen::GeneralOracle, rm::Model, prefs::Dict{Symbol,Any})

    # Extract preferences we care about
    gen.use_cuts          = get(prefs, :prefer_cuts, false)
    gen.cut_tol           = get(prefs, :cut_tol, 1e-6)
    gen.debug_printcut    = get(prefs, :debug_printcut, false)

    rd = getRobust(rm)

    #-------------------------------------------------------------------
    # Cutting plane setup
    if gen.use_cuts || prefs[:active_cuts]
        # Create an LP/SOCP that we'll use to solve the cut problem
        # Copy the uncertainty set from the original problem
        gen.cut_model = Model()
        gen.cut_model.solver   = isa(rd.cutsolver,JuMP.UnsetSolver) ? rm.solver : rd.cutsolver
        gen.cut_model.numCols  = rd.numUncs
        gen.cut_model.colNames = rd.uncNames
        gen.cut_model.colLower = rd.uncLower
        gen.cut_model.colUpper = rd.uncUpper
        gen.cut_model.colCat   = rd.uncCat
        gen.cut_vars = [Variable(gen.cut_model, i) for i in 1:rd.numUncs]

        # Covert polyhedral constraints
        for c in rd.uncertaintyset
            push!(gen.cut_model.linconstr, 
                LinearConstraint(uaff_to_aff(c.terms,gen.cut_vars), c.lb, c.ub))
        end
        # Norm constraints
        for norm_c in rd.normconstraints
            if isa(norm_c, UncNormConstraint{2})
                # Ellispoidal constraint
                # Input: ‖[a₁ᵀu, a₂ᵀu, ...]‖₂ ≤ Γ
                # Output: yᵢ = aᵢᵀu, t = Γ
                #         Σyᵢ^2 ≤ t^2
                normexp = norm_c.normexpr
                terms   = normexp.norm.terms
                rhs     = -normexp.aff.constant / normexp.coeff
                n_terms = length(terms)
                @defVar(gen.cut_model, y[1:n_terms])
                for i in 1:n_terms
                    @addConstraint(gen.cut_model, y[i] ==
                        uaff_to_aff(terms[i],gen.cut_vars))
                end
                @defVar(gen.cut_model, t == rhs)
                @addConstraint(gen.cut_model, dot(y,y) <= t^2)
            elseif isa(norm_c, UncNormConstraint{1})
                # L1 norm constraint
                # Input: ‖[a₁ᵀu, a₂ᵀu, ...]‖₁ ≤ Γ
                # Output: yᵢ ≥ aᵢᵀu, yᵢ ≥ -aᵢᵀu
                #         ∑yᵢ ≤ Γ
                normexp = norm_c.normexpr
                terms   = normexp.norm.terms
                rhs     = -normexp.aff.constant / normexp.coeff
                n_terms = length(terms)
                @defVar(gen.cut_model, y[1:n_terms])
                for i in 1:n_terms
                    @addConstraint(gen.cut_model, y[i] ≥
                         uaff_to_aff(terms[i],gen.cut_vars))
                    @addConstraint(gen.cut_model, y[i] ≥
                        -uaff_to_aff(terms[i],gen.cut_vars))
                end
                @addConstraint(gen.cut_model, sum(y) <= rhs)
            elseif isa(norm_c, UncNormConstraint{Inf})
                # L∞ norm constraint
                # Input: ‖[a₁ᵀu, a₂ᵀu, ...]‖∞ ≤ Γ
                # Output: aᵢᵀu ≤ Γ, aᵢᵀu ≥ -Γ
                normexp = norm_c.normexpr
                terms   = normexp.norm.terms
                rhs     = -normexp.aff.constant / normexp.coeff
                n_terms = length(terms)
                for i in 1:n_terms
                    @addConstraint(gen.cut_model, 
                        uaff_to_aff(terms[i],gen.cut_vars) ≤ +rhs)
                    @addConstraint(gen.cut_model,
                        uaff_to_aff(terms[i],gen.cut_vars) ≥ -rhs)
                end
            else
                error("Unrecognized norm in uncertainty set!")
            end
        end
    end

    #-------------------------------------------------------------------
    # Reformulation setup
    if !gen.use_cuts
        # Statement of primal-dual pair:
        # PRIMAL
        # max  cᵀx
        #  st  A x     ≤ b   ::  π       [  #(linear constraints) ]
        #      ‖Fx+f‖₂ ≤ Γ₂  ::  β,β′    [ 1 + #(terms in 2-norm) ]
        #      ‖Gx+g‖₁ ≤ Γ₁  ::  α,α′    [ 1 + #(terms in 1-norm) ]
        #      ‖Hx+h‖∞ ≤ Γ∞  ::  ω,ω′    [ 2 × #(terms in ∞-norm) ]
        #
        # DUAL
        #          {  ELLIPSE } {    L1 NORM     } {      L∞ NORM      }
        # min  bᵀπ + fᵀβ + Γ₂β′ + gᵀα + (gᵀ1+Γ₁)α′ + (Γ∞-hᵀω) - (Γ∞+hᵀω)
        #  st  Aᵀπ - Fᵀβ        - Gᵀα - (Gᵀ1   )α′ +     Hᵀω  +     Hᵀω′ = c
        #      π  ≥ 0
        #      β′ ≥ ‖β‖₂  ('ell_lhs' ≥ 'ell_rhs')
        #      α′ ≥ -½α   ( 'l1_lhs' ≥  'l1_rhs')
        #      α ≤ 0, α′ ≥ 0
        #      ω ≥ 0, ω′ ≤ 0
        #
        # In this setup phase we build the structure of the dual but do not
        # attach it to the model. Rather, for each constraint we will spawn
        # a new set of variables and constraints according to the structure
        # we determine here.
        
        # Number of dual variables
        # =   Number of linear constraints 
        #  +  Number of terms in 2-norm constraints + 1
        #  +  Number of terms in 1-norm constraints + 1
        #  +  Number of terms in ∞-norm constraints × 2
        #  +  Number of bounds on uncertain parameters (later)
        num_dualvar = length(rd.uncertaintyset)
        for norm_c in rd.normconstraints
            if isa(norm_c, UncNormConstraint{2}) || isa(norm_c, UncNormConstraint{1})
                num_dualvar += length(norm_c.normexpr.norm.terms) + 1
            elseif isa(norm_c, UncNormConstraint{Inf})
                num_dualvar += length(norm_c.normexpr.norm.terms) * 2
            else
                error("Unrecognized norm in uncertainty set!")
            end
        end
        # Number of linear dual constraints = number of uncertain parameters
        num_dualcon  = rd.numUncs
        # Store Aᵀ as row-wise sparse vectors
        dual_A       = [(@compat Tuple{Int,Float64})[] for i in 1:rd.numUncs]
        # Store dual objective coefficients
        dual_objs    = zeros(num_dualvar)
        # Store dual variable sense
        dual_vartype = fill(:(>=), num_dualvar)
        # Store dual linear constraint sense
        # All are at equality because all bounds on the uncertain parameters
        # are converted to constraints
        dual_contype = fill(:(==), num_dualcon)

        # Linear constraints form Aᵀ, and set vartype and objective of π
        for (aff_ind, aff_con) in enumerate(rd.uncertaintyset)
            lhs = aff_con.terms
            for (coef,uncp) in lhs
                push!(dual_A[uncp.id], (aff_ind, coef))
            end
            dual_objs[aff_ind]    =   rhs(aff_con)
            dual_vartype[aff_ind] = sense(aff_con)
        end

        # Set ellipse objective coefficients for each β.β′
        # Track index of last set dual - in this case dual for the
        # last linear constraint. We'll update after each ellipse.
        start_ind = length(rd.uncertaintyset)
        for norm_c in rd.normconstraints
            if isa(norm_c, UncNormConstraint{2})
                # Extract fields from norm constraint
                normexp = norm_c.normexpr
                terms   = normexp.norm.terms
                rhs     = -normexp.aff.constant / normexp.coeff

                # Dual β, the RHS of β′ ≥ ‖β‖
                dual_ell_rhs_idx = Int[]
                for (term_ind, term) in enumerate(terms)
                    dual_ind = start_ind + term_ind
                    dual_objs[dual_ind]    = term.constant
                    dual_vartype[dual_ind] = :Qrhs
                    # Track indices for this ellipse
                    push!(dual_ell_rhs_idx, dual_ind)
                end
                push!(gen.dual_ell_rhs_idxs, dual_ell_rhs_idx)

                # Dual β′, the LHS of β′ ≥ ‖β‖
                dual_ind = start_ind + length(terms) + 1
                dual_objs[dual_ind]    = rhs
                dual_vartype[dual_ind] = :Qlhs
                push!(gen.dual_ell_lhs_idxs, dual_ind)
                start_ind += length(terms) + 1
            end
        end

        # Same as above, but for 1-norm constraints
        for norm_c in rd.normconstraints
            if isa(norm_c, UncNormConstraint{1})
                # Extract fields from norm constraint
                normexp = norm_c.normexpr
                terms   = normexp.norm.terms
                rhs     = -normexp.aff.constant / normexp.coeff

                # Dual α, the RHS of α′ ≥ -½α
                total_g = 0.0
                dual_l1_rhs_idx = Int[]
                for (term_ind, term) in enumerate(terms)
                    dual_ind = start_ind + term_ind
                    total_g += term.constant
                    dual_objs[dual_ind]    = term.constant
                    dual_vartype[dual_ind] = :L1rhs
                    # Track indices for this ellipse
                    push!(dual_l1_rhs_idx, dual_ind)
                end
                push!(gen.dual_l1_rhs_idxs, dual_l1_rhs_idx)

                # Dual β′, the LHS of β′ ≥ ‖β‖
                dual_ind = start_ind + length(terms) + 1
                dual_objs[dual_ind]    = total_g + rhs
                dual_vartype[dual_ind] = :L1lhs
                push!(gen.dual_l1_lhs_idxs, dual_ind)
                start_ind += length(terms) + 1
            end
        end

        # Same as above, but for ∞-norm constraints
        for norm_c in rd.normconstraints
            if isa(norm_c, UncNormConstraint{Inf})
                # Extract fields from norm constraint
                normexp = norm_c.normexpr
                terms   = normexp.norm.terms
                rhs     = -normexp.aff.constant / normexp.coeff

                # Dual ω and ω′
                dual_ω_idx, dual_ω′_idx = Int[], Int[]
                for (term_ind, term) in enumerate(terms)
                    # First do ωᵢ
                    dual_ind = start_ind + term_ind
                    dual_objs[dual_ind]    =+rhs - term.constant
                    dual_vartype[dual_ind] = :(ω)
                    push!(dual_ω_idx, dual_ind)
                    # Then do ω′ᵢ
                    dual_ind = start_ind + term_ind + length(terms)
                    dual_objs[dual_ind]    = -rhs - term.constant
                    dual_vartype[dual_ind] = :(ω′)
                    push!(dual_ω′_idx, dual_ind)
                end
                push!(gen.dual_ω_idxs,  dual_ω_idx)
                push!(gen.dual_ω′_idxs, dual_ω′_idx)
                start_ind += length(terms) * 2
            end
        end

        # Bounds on uncertain parameters are handled as linear constraints
        for i in 1:rd.numUncs
            L, U = rd.uncLower[i], rd.uncUpper[i]
            # If it has a lower bound...
            if L != -Inf
                num_dualvar += 1
                push!(dual_A[i],    (num_dualvar, 1.0))
                push!(dual_objs,    L)
                push!(dual_vartype, :(>=))
            end
            # If it has an upper bound...
            if U != +Inf
                num_dualvar += 1
                push!(dual_A[i],    (num_dualvar, 1.0))
                push!(dual_objs,    U)
                push!(dual_vartype, :(<=))
            end
        end

        # Store the reformulation in the oracle
        gen.num_dualvar  = num_dualvar
        gen.dual_A       = dual_A
        gen.dual_objs    = dual_objs
        gen.dual_vartype = dual_vartype
        gen.dual_contype = dual_contype
    end  # end reformulation preparation
end
Base.precompile(setup, (GeneralOracle, Model, Dict{Symbol,Any}))

#-----------------------------------------------------------------------
# Cutting planes
# Given constraints with uncertain parameters, create new constraints
# if they would cause the current solution to be infeasible

function generateCut(gen::GeneralOracle, master::Model, rm::Model, inds::Vector{Int}, active=false)
    # If not doing cuts...
    (!gen.use_cuts && !active) && return Any[]

    rd = getRobust(rm)
    master_sol = master.colVal
    new_cons = Any[]

    for con_ind in inds
        con = get_uncertain_constraint(rm, con_ind)

        # Update the cutting plane problem's objective, and solve
        cut_sense, unc_obj_coeffs, lhs_const = JuMPeR.build_cut_objective_sparse(con, master_sol)
        @setObjective(gen.cut_model, cut_sense, sum{u[2]*gen.cut_vars[u[1]], u=unc_obj_coeffs})
        cut_solve_status = solve(gen.cut_model, suppress_warnings=true)
        cut_solve_status != :Optimal &&
            error("GeneralOracle: cutting plane problem is infeasible or unbounded!")
        lhs_of_cut = getObjectiveValue(gen.cut_model) + lhs_const

        # SUBJECT TO CHANGE: active cut detection
        if active
            push!(rd.activecuts[con_ind], 
                cut_to_scen(gen.cut_model.colVal, 
                    check_cut_status(con, lhs_of_cut, gen.cut_tol) == :Active))
            continue
        end
        
        # Check violation
        if check_cut_status(con, lhs_of_cut, gen.cut_tol) != :Violate
            gen.debug_printcut && debug_printcut(rm,master,gen,lhs_of_cut,con,nothing)
            continue  # No violation, no new cut
        end
        
        # Create and add the new constraint
        new_con = JuMPeR.build_certain_constraint(master, con, gen.cut_model.colVal)
        gen.debug_printcut && debug_printcut(rm,master,gen,lhs_of_cut,con,new_con)
        push!(new_cons, new_con)
    end
    
    return new_cons
end
Base.precompile(generateCut, (GeneralOracle, Model, Model, Vector{Int}, Bool))


function debug_printcut(rm,m,w,lhs,con,new_con)
    println("BEGIN DEBUG :debug_printcut")
    #convert_model!(con, rm)
    println("  Constraint:  ", con)
    #convert_model!(con, m)
    println("  Master sol:  ")
    for j in 1:length(m.colNames)
        println("    ", m.colNames[j], "  ", m.colVal[j])
    end
    print(w.cut_model)
    #println("  Solve status ", cut_solve_status)
    println("  Cut sol:     ")
    for j in 1:length(w.cut_model.colNames)
        println("    ", w.cut_model.colNames[j], "  ", w.cut_model.colVal[j])
    end
    println("  OrigLHS val: ", lhs)
    println("  Sense:       ", sense(con))
    println("  con.lb/ub:   ", con.lb, "  ", con.ub)
    println("  new con:  ", new_con)
    println("END DEBUG   :debug_printcut")
end




#-----------------------------------------------------------------------
# Reformulation
# In setup() the 'template' for reformulations was calculated. We now
# apply this template to every constraint that the user would like to
# reformulate.

# generateReform
# Modifies the master problem in place, returning the number of
# constraints that were reformulated
function generateReform(gen::GeneralOracle, master::Model, rm::Model, inds::Vector{Int})
    # If not doing reform...
    gen.use_cuts && return 0
    # Apply the reformulation to all relevant constraints
    for ind in inds
        apply_reform(gen, master, rm, ind)
    end
    return length(inds)
end
Base.precompile(generateReform, (GeneralOracle, Model, Model, Vector{Int}))

# apply_reform
# Does the hard work for a given constraint
function apply_reform(gen::GeneralOracle, master::Model, rm::Model, con_ind::Int)
    rd = getRobust(rm)
    con = get_uncertain_constraint(rm, con_ind)

    # Pull the 'template' out of the oracle for easier access
    num_dualvar  = gen.num_dualvar
    num_dualcon  = rd.numUncs
    dual_A       = gen.dual_A
    dual_objs    = gen.dual_objs
    dual_vartype = gen.dual_vartype
    dual_contype = gen.dual_contype 

    # Initialize the affine expressions that are equal to the coefficients
    # of the uncertain parameters in the primal of the cutting plane
    # problem, and RHS values for the dual problem
    dual_rhs    = [AffExpr() for i in 1:num_dualcon]
    # This is constraint that will replace the existing uncertain constraint
    new_lhs     = AffExpr()
    # We do all reformulation as if the constraint is a <= constraint
    # This necessitates mulitplying through by -1 if it is a >= constraint
    sign_flip   = sense(con) == :(<=) ? +1.0 : -1.0
    # Initialize the RHS of the new constraint
    new_rhs     = rhs(con) * sign_flip
    # Extract the LHS of the constraint
    orig_lhs    = con.terms

    # Original LHS terms are of form (aᵢᵀu + bᵢ) xᵢ. 
    # Collect the certain terms (bᵢxᵢ) of the uncertain constraint,
    # and append them to the new LHS directly
    for (uaff,var) in orig_lhs
        if uaff.constant != 0.0
            push!(new_lhs, uaff.constant * sign_flip,
                    Variable(master, var.col))
        end
    end
    
    # Rearrange from ∑ᵢ (aᵢᵀu) xᵢ to ∑ⱼ (cⱼᵀx) uⱼ, as the cⱼᵀx 
    # are the RHS of the dual. While constructing, we check for
    # integer uncertain parameters, which we cannot reformulate
    for (uaff,var) in orig_lhs
        for (coeff, uncparam) in uaff
            getCategory(uncparam) != :Cont &&
                error("Integer uncertain parameters not supported in reformulation.")
            push!(dual_rhs[uncparam.id], coeff * sign_flip,
                    Variable(master, var.col))
        end
    end
    # We also need the standalone aᵀu not related to any variable
        for (coeff, uncparam) in orig_lhs.constant
            getCategory(uncparam) != :Cont &&
                error("Integer uncertain parameters not supported in reformulation.")
            dual_rhs[uncparam.id].constant += coeff * sign_flip
        end

    # Create dual variables for this specific contraint and add
    # them to the new constraint's LHS
    dual_vars = Variable[]
    for ind in 1:num_dualvar
        # Free:   equality,     RHS of ellipse
        # NonNeg: less-than,    LHS of ellipse, LHS of L1 norm, ω
        # NonPos: greater-than,                 RHS of L1 norm, ω′
        vt = dual_vartype[ind]
        lbound = (vt == :(<=) || vt == :Qlhs || vt == :L1lhs || vt == :ω ) ? 0 : -Inf
        ubound = (vt == :(>=) ||                vt == :L1rhs || vt == :ω′) ? 0 : +Inf
        vname = "π"
        vt == :Qlhs  && (vname="β′")
        vt == :Qrhs  && (vname="β" )
        vt == :L1lhs && (vname="α′")
        vt == :L1rhs && (vname="α" )
        vt == :ω     && (vname="ω")
        vt == :ω′    && (vname="ω′" )
        new_v = Variable(master,lbound,ubound,:Cont,"_$(vname)_$(con_ind)_$(ind)")
        push!(dual_vars, new_v)
        push!(new_lhs, dual_objs[ind], new_v)
    end

    # Add the new deterministic constraint to the master problem
    @addConstraint(master, new_lhs <= new_rhs)

    # Add the additional new constraints
    for unc in 1:rd.numUncs
        new_lhs = AffExpr()
        # aᵀπ
        for (ind,coeff) in dual_A[unc]
            push!(new_lhs, coeff, dual_vars[ind])
        end

        # Norms
        ell_idx, l1_idx, l∞_idx = 0, 0, 0
        for norm_c in rd.normconstraints
            # 2-norm    -Fᵀβ
            if isa(norm_c, UncNormConstraint{2})
                ell_idx += 1
                terms = norm_c.normexpr.norm.terms
                for (term_ind, term) in enumerate(terms)
                    for (coeff,uncparam) in term
                        # Is it a match?
                        uncparam.id != unc && continue
                        # F ≠ 0 for this uncertain parameter and term
                        ell_rhs_idxs = gen.dual_ell_rhs_idxs[ell_idx]
                        push!(new_lhs, -coeff, dual_vars[ell_rhs_idxs[term_ind]])
                    end
                end
            # 1-norm    -Gᵀα -(Gᵀ1)α′
            elseif isa(norm_c, UncNormConstraint{1})
                l1_idx += 1
                terms = norm_c.normexpr.norm.terms
                for (term_ind, term) in enumerate(terms)
                    for (coeff,uncparam) in term
                        # Is it a match?
                        uncparam.id != unc && continue
                        # G ≠ 0 for this uncertain parameter and term
                        push!(new_lhs, -coeff, dual_vars[gen.dual_l1_rhs_idxs[l1_idx][term_ind]])
                        push!(new_lhs, -coeff, dual_vars[gen.dual_l1_lhs_idxs[l1_idx]          ])
                    end
                end
            # ∞-norm    Hᵀω + Hᵀω′
            elseif isa(norm_c, UncNormConstraint{Inf})
                l∞_idx += 1
                terms = norm_c.normexpr.norm.terms
                for (term_ind, term) in enumerate(terms)
                    for (coeff,uncparam) in term
                        # Is it a match?
                        uncparam.id != unc && continue
                        # H ≠ 0 for this uncertain parameter and term
                        push!(new_lhs, coeff, dual_vars[gen.dual_ω_idxs[ l∞_idx][term_ind]])
                        push!(new_lhs, coeff, dual_vars[gen.dual_ω′_idxs[l∞_idx][term_ind]])
                    end
                end
            end
        end

        ct = dual_contype[unc]
        ct == :(==) && @addConstraint(master, new_lhs == dual_rhs[unc])
        ct == :(<=) && @addConstraint(master, new_lhs <= dual_rhs[unc])
        ct == :(>=) && @addConstraint(master, new_lhs >= dual_rhs[unc])
    end


    ell_idx,  l1_idx = 0, 0
    for norm_c in rd.normconstraints
        # Impose β′ ≥ ‖β‖₂
        if isa(norm_c, UncNormConstraint{2})
            ell_idx += 1
            β′ = dual_vars[gen.dual_ell_lhs_idxs[ell_idx]]
            β  = dual_vars[gen.dual_ell_rhs_idxs[ell_idx]]
            @addConstraint(master, dot(β,β) <= β′*β′)
        # Impose α′ ≥ -½α
        elseif isa(norm_c, UncNormConstraint{1})
            l1_idx += 1
            α′ = dual_vars[gen.dual_l1_lhs_idxs[l1_idx]]
            α  = dual_vars[gen.dual_l1_rhs_idxs[l1_idx]]
            for αᵢ in α
                @addConstraint(master, α′ >= -0.5αᵢ)
            end
        end
    end

    return true
end
Base.precompile(apply_reform, (GeneralOracle, Model, Model, Int))