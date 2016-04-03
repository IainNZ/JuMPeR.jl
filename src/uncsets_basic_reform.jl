#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2016: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# src/uncsets_basic_reform.jl
# Reformulation-specific functionality for the BasicUncertaintySet,
# including setup_set_reform and generate_reform.
# Included by src/uncsets_basic.jl
#-----------------------------------------------------------------------


"""
    setup_sets_reform(BasicUncertaintySet, RobustModel)

Prepares structures for reformulating with the BasicUncertaintySet.
"""
function setup_set_reform(us::BasicUncertaintySet, rm::Model)
    # Statement of primal-dual pair:
    # PRIMAL
    # max  cᵀx
    #  st  A x     ≤ b   ::  π       [  #(linear constraints) ]
    #      ‖Fx+f‖₁ ≤ Γ₁  ::  α,α′    [ 1 + #(terms in 1-norm) ]
    #      ‖Gx+g‖₂ ≤ Γ₂  ::  β,β′    [ 1 + #(terms in 2-norm) ]
    #      ‖Hx+h‖∞ ≤ Γ∞  ::  ω,ω′    [ 2 × #(terms in ∞-norm) ]
    #
    # DUAL
    #           {    L1 NORM     } {  ELLIPSE } {      L∞ NORM      }
    # min  bᵀπ + fᵀα + (fᵀ1+Γ₁)α′ + gᵀβ + Γ₂β′ + (Γ∞-hᵀω) - (Γ∞+hᵀω)
    #  st  Aᵀπ - Fᵀα - (Fᵀ1   )α′ - Gᵀβ        +     Hᵀω  +     Hᵀω′ = c
    #      π  ≥ 0
    #      α′ ≥ -½α   ( 'l1_lhs' ≥  'l1_rhs')
    #      β′ ≥ ‖β‖₂  ('ell_lhs' ≥ 'ell_rhs')
    #      α ≤ 0, α′ ≥ 0
    #      ω ≥ 0, ω′ ≤ 0
    #
    # In this setup phase we build the structure of the dual but do not
    # attach it to the model. Rather, for each constraint we will spawn
    # a new set of variables and constraints according to the structure
    # we determine here.

    # Number of dual variables
    # =   Number of linear constraints
    #  +  Number of terms in 1-norm constraints + 1
    #  +  Number of terms in 2-norm constraints + 1
    #  +  Number of terms in ∞-norm constraints × 2
    #  +  Number of bounds on uncertain parameters (later)
    rmext = get_robust(rm)
    num_dualvar = length(us.linear_constraints)
    for norm_c in us.norm_constraints
        if isa(norm_c, UncSetNormConstraint{2}) || isa(norm_c, UncSetNormConstraint{1})
            num_dualvar += length(norm_c.normexpr.norm.terms) + 1
        elseif isa(norm_c, UncSetNormConstraint{Inf})
            num_dualvar += length(norm_c.normexpr.norm.terms) * 2
        else
            error("Unrecognized norm in uncertainty set!")
        end
    end
    # Number of linear dual constraints = number of uncertain parameters
    num_dualcon  = rmext.num_uncs
    # Store Aᵀ as row-wise sparse vectors
    dual_A       = [Tuple{Int,Float64}[] for i in 1:rmext.num_uncs]
    # Store dual objective coefficients
    dual_objs    = zeros(num_dualvar)
    # Store dual variable sense
    dual_vartype = fill(:(>=), num_dualvar)
    # Store dual linear constraint sense
    # All are at equality because all bounds on the uncertain parameters
    # are converted to constraints
    dual_contype = fill(:(==), num_dualcon)

    # Linear constraints form Aᵀ, and set vartype and objective of π
    for (aff_ind, aff_con) in enumerate(us.linear_constraints)
        lhs = aff_con.terms
        for (coef,uncp) in linearterms(lhs)
            push!(dual_A[uncp.id], (aff_ind, coef))
        end
        dual_objs[aff_ind]    =   rhs(aff_con)
        dual_vartype[aff_ind] = sense(aff_con)
    end

    # Set ellipse objective coefficients for each β.β′
    # Track index of last set dual - in this case dual for the
    # last linear constraint. We'll update after each ellipse.
    start_ind = length(us.linear_constraints)
    for norm_c in us.norm_constraints
        if isa(norm_c, UncSetNormConstraint{2})
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
            push!(us.dual_ell_rhs_idxs, dual_ell_rhs_idx)

            # Dual β′, the LHS of β′ ≥ ‖β‖
            dual_ind = start_ind + length(terms) + 1
            dual_objs[dual_ind]    = rhs
            dual_vartype[dual_ind] = :Qlhs
            push!(us.dual_ell_lhs_idxs, dual_ind)
            start_ind += length(terms) + 1
        end
    end

    # Same as above, but for 1-norm constraints
    for norm_c in us.norm_constraints
        if isa(norm_c, UncSetNormConstraint{1})
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
            push!(us.dual_l1_rhs_idxs, dual_l1_rhs_idx)

            # Dual β′, the LHS of β′ ≥ ‖β‖
            dual_ind = start_ind + length(terms) + 1
            dual_objs[dual_ind]    = total_g + rhs
            dual_vartype[dual_ind] = :L1lhs
            push!(us.dual_l1_lhs_idxs, dual_ind)
            start_ind += length(terms) + 1
        end
    end

    # Same as above, but for ∞-norm constraints
    for norm_c in us.norm_constraints
        if isa(norm_c, UncSetNormConstraint{Inf})
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
            push!(us.dual_ω_idxs,  dual_ω_idx)
            push!(us.dual_ω′_idxs, dual_ω′_idx)
            start_ind += length(terms) * 2
        end
    end

    # Bounds on uncertain parameters are handled as linear constraints
    for i in 1:rmext.num_uncs
        L, U = rmext.unc_lower[i], rmext.unc_upper[i]
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

    # Store the reformulation
    us.num_dualvar  = num_dualvar
    us.dual_A       = dual_A
    us.dual_objs    = dual_objs
    us.dual_vartype = dual_vartype
    us.dual_contype = dual_contype
end


"""
    generate_reform(BasicUncertaintySet, ...)

Modifies the problem in place.
"""
function generate_reform(us::BasicUncertaintySet, rm::Model, idxs::Vector{Int})
    # If not doing reform...
    if us.use_cuts
        return 0
    end
    # Apply the reformulation to all relevant constraints
    for idx in idxs
        apply_reform(us, rm, idx)
    end
end


"""
    apply_reform(BasicUncertaintySet, ...)

Reformulates a single constraint.
"""
function apply_reform(us::BasicUncertaintySet, rm::Model, idx::Int)
    rmext = get_robust(rm)::RobustModelExt
    con = rmext.unc_constraints[idx]

    # Pull the 'template' out for easier access
    num_dualvar  = us.num_dualvar
    num_dualcon  = rmext.num_uncs
    dual_A       = us.dual_A
    dual_objs    = us.dual_objs
    dual_vartype = us.dual_vartype
    dual_contype = us.dual_contype

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
    for (uaff,var) in linearterms(orig_lhs)
        if uaff.constant != 0.0
            push!(new_lhs, uaff.constant * sign_flip, var)
        end
    end

    # Rearrange from ∑ᵢ (aᵢᵀu) xᵢ to ∑ⱼ (cⱼᵀx) uⱼ, as the cⱼᵀx
    # are the RHS of the dual. While constructing, we check for
    # integer uncertain parameters, which we cannot reformulate
    for (uaff,var) in linearterms(orig_lhs)
        for (coeff, uncparam) in linearterms(uaff)
            rmext.unc_cat[uncparam.id] != :Cont &&
                error("Integer uncertain parameters not supported in reformulation.")
            push!(dual_rhs[uncparam.id], coeff * sign_flip, var)
        end
    end
    # We also need the standalone aᵀu not related to any variable
        for (coeff, uncparam) in linearterms(orig_lhs.constant)
            rmext.unc_cat[uncparam.id] != :Cont &&
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
        new_v = Variable(rm,lbound,ubound,:Cont,"_$(vname)_$(idx)_$(ind)")
        push!(dual_vars, new_v)
        push!(new_lhs, dual_objs[ind], new_v)
    end

    # Add the new deterministic constraint to the problem
    @addConstraint(rm, new_lhs <= new_rhs)

    # Add the additional new constraints
    for unc in 1:rmext.num_uncs
        new_lhs = AffExpr()
        # aᵀπ
        for (ind,coeff) in dual_A[unc]
            push!(new_lhs, coeff, dual_vars[ind])
        end

        # Norms
        ell_idx, l1_idx, l∞_idx = 0, 0, 0
        for norm_c in us.norm_constraints
            # 2-norm    -Fᵀβ
            if isa(norm_c, UncSetNormConstraint{2})
                ell_idx += 1
                terms = norm_c.normexpr.norm.terms
                for (term_ind, term) in enumerate(terms)
                    for (coeff,uncparam) in linearterms(term)
                        # Is it a match?
                        uncparam.id != unc && continue
                        # F ≠ 0 for this uncertain parameter and term
                        ell_rhs_idxs = us.dual_ell_rhs_idxs[ell_idx]
                        push!(new_lhs, -coeff, dual_vars[ell_rhs_idxs[term_ind]])
                    end
                end
            # 1-norm    -Gᵀα -(Gᵀ1)α′
            elseif isa(norm_c, UncSetNormConstraint{1})
                l1_idx += 1
                terms = norm_c.normexpr.norm.terms
                for (term_ind, term) in enumerate(terms)
                    for (coeff,uncparam) in linearterms(term)
                        # Is it a match?
                        uncparam.id != unc && continue
                        # G ≠ 0 for this uncertain parameter and term
                        push!(new_lhs, -coeff, dual_vars[us.dual_l1_rhs_idxs[l1_idx][term_ind]])
                        push!(new_lhs, -coeff, dual_vars[us.dual_l1_lhs_idxs[l1_idx]          ])
                    end
                end
            # ∞-norm    Hᵀω + Hᵀω′
            elseif isa(norm_c, UncSetNormConstraint{Inf})
                l∞_idx += 1
                terms = norm_c.normexpr.norm.terms
                for (term_ind, term) in enumerate(terms)
                    for (coeff,uncparam) in linearterms(term)
                        # Is it a match?
                        uncparam.id != unc && continue
                        # H ≠ 0 for this uncertain parameter and term
                        push!(new_lhs, coeff, dual_vars[us.dual_ω_idxs[ l∞_idx][term_ind]])
                        push!(new_lhs, coeff, dual_vars[us.dual_ω′_idxs[l∞_idx][term_ind]])
                    end
                end
            end
        end

        ct = dual_contype[unc]
        ct == :(==) && @addConstraint(rm, new_lhs == dual_rhs[unc])
        ct == :(<=) && @addConstraint(rm, new_lhs <= dual_rhs[unc])
        ct == :(>=) && @addConstraint(rm, new_lhs >= dual_rhs[unc])
    end


    ell_idx,  l1_idx = 0, 0
    for norm_c in us.norm_constraints
        # Impose β′ ≥ ‖β‖₂
        if isa(norm_c, UncSetNormConstraint{2})
            ell_idx += 1
            β′ = dual_vars[us.dual_ell_lhs_idxs[ell_idx]]
            β  = dual_vars[us.dual_ell_rhs_idxs[ell_idx]]
            @addConstraint(rm, dot(β,β) <= β′*β′)
        # Impose α′ ≥ -½α
        elseif isa(norm_c, UncSetNormConstraint{1})
            l1_idx += 1
            α′ = dual_vars[us.dual_l1_lhs_idxs[l1_idx]]
            α  = dual_vars[us.dual_l1_rhs_idxs[l1_idx]]
            for αᵢ in α
                @addConstraint(rm, α′ >= -0.5αᵢ)
            end
        end
    end

    return true
end
