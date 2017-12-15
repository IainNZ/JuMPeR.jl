#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2016: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# src/print.jl
# All "pretty printers" for JuMPeR types.
#-----------------------------------------------------------------------

import JuMP: REPLMode, IJuliaMode, PrintSymbols
import JuMP: repl, ijulia
import JuMP: PRINT_ZERO_TOL, DIMS
import JuMP: str_round, getmeta
import JuMP: aff_str, aff_str
import JuMP: cont_str, con_str

# helper to look up corresponding JuMPContainerData
printdata(v::JuMPContainer{Uncertain}) = get_robust(getmeta(v, :model)).uncData[v]
function printdata(u::Array{Uncertain})
    if isempty(u)
        error("Cannot locate printing data for an empty array")
    end
    m = first(u).m
    get_robust(m).uncData[u]
end


#------------------------------------------------------------------------
## RobustModel
#------------------------------------------------------------------------
# Called by the JuMP print hook
function print_robust(io::IO, m::Model)
    rmext = get_robust(m)
    # First, display normal model stuff
    JuMP.print(io, m, ignore_print_hook=true)
    # Now print adaptive variable info
    for i in 1:rmext.num_adps
        adp_name = adp_str(REPLMode,m,i) * "(ξ)"
        adp_lb, adp_ub = rmext.adp_lower[i], rmext.adp_upper[i]
        str_lb, str_ub = str_round(adp_lb), str_round(adp_ub)
        adp_cat = rmext.adp_cat[i]
        if adp_cat == :Bin  # x binary
            str = string(adp_name, " ", repl[:in], " ",
                        repl[:open_set], "0,1", repl[:close_set])
        elseif adp_lb == -Inf && adp_ub == +Inf # Free
            str = string(adp_name, " free")
        elseif adp_lb == -Inf  # No lower bound
            str = string(adp_name, " ", repl[:leq], " ", str_ub)
        elseif adp_ub == +Inf  # No upper bound
            str = string(adp_name, " ", repl[:geq], " ", str_lb)
        else
            str = string(str_lb, " ", repl[:leq], " ", adp_name,
                            " ", repl[:leq], " ", str_ub)
        end
        if adp_cat == :Int
            str *= string(", ", repl[:integer])
        end
        println(io, " ", str)
    end

    # Now print stuff particular to JuMPeR
    println(io, "Uncertain constraints:")
    for uc in rmext.unc_constraints
        println(io, uc)
    end
    #println(io, "Uncertainty set:")
    #for us in rmext.uncertaintyset
    #    println(io, us)
    #end
    #for nc in rmext.normconstraints
    #    println(io, nc)
    #end
    println(io, "Uncertain parameters:")
    # Display indexed uncertains
    in_dictlist = falses(rmext.num_uncs)
    for d in rmext.dictList
        println(io, cont_str(REPLMode,d))
        for it in JuMP._values(d)  # Mark uncertains in JuMPContainer as printed
            in_dictlist[it.id] = true
        end
    end

    # Display non-indexed variables
    for i in 1:rmext.num_uncs
        in_dictlist[i] && continue
        str = ""
        unc_name = unc_str(REPLMode,m,i)
        unc_lb, unc_ub = rmext.unc_lower[i], rmext.unc_upper[i]
        str_lb, str_ub = str_round(unc_lb), str_round(unc_ub)
        unc_cat = rmext.unc_cat[i]
        if unc_cat == :Bin  # x binary
            str = string(unc_name, " ", repl[:in], " ",
                        repl[:open_set], "0,1", repl[:close_set])
        elseif unc_lb == -Inf && unc_ub == +Inf # Free
            str = string(unc_name, " free")
        elseif unc_lb == -Inf  # No lower bound
            str = string(unc_name, " ", repl[:leq], " ", str_ub)
        elseif unc_ub == +Inf  # No upper bound
            str = string(unc_name, " ", repl[:geq], " ", str_lb)
        else
            str = string(str_lb, " ", repl[:leq], " ", unc_name,
                            " ", repl[:leq], " ", str_ub)
        end
        if unc_cat == :Int
            str *= string(", ", repl[:integer])
        end
        println(io, str)
    end
end

#------------------------------------------------------------------------
## Uncertain
#------------------------------------------------------------------------
Base.print(io::IO, u::Uncertain) = print(io, unc_str(REPLMode,u))
Base.show( io::IO, u::Uncertain) = print(io, unc_str(REPLMode,u))
#Base.writemime(io::IO, ::MIME"text/latex", u::Uncertain) =
#    print(io, unc_str(IJuliaMode,u,mathmode=false))
function unc_str(mode, m::Model, id::Int)
    rmext = get_robust(m)
    uncNames = rmext.unc_names
    if uncNames[id] == ""
        for cont in rmext.dictList
            fill_unc_names(mode, uncNames, cont)
        end
    end
    return uncNames[id] == "" ? "unc_$id" : uncNames[id]
end
function fill_unc_names{N}(mode, uncNames, u::JuMP.JuMPArray{Uncertain,N})
    data = printdata(u)
    idxsets = data.indexsets
    lengths = map(length, idxsets)
    name = data.name
    cprod = cumprod([lengths...])
    for (ind,unc) in enumerate(u.innerArray)
        idx_strs = [string( idxsets[1][mod1(ind,lengths[1])] )]
        for i = 2:N
            push!(idx_strs, string(idxsets[i][Int(ceil(mod1(ind,cprod[i]) / cprod[i-1]))]))
        end
        #if mode == IJuliaMode
        #    uncNames[unc.id] = string(name, "_{", join(idx_strs,",") , "}")
        #else
            uncNames[unc.id] = string(name,  "[", join(idx_strs,",") , "]")
        #end
    end
end
function fill_unc_names(mode, uncNames, u::JuMP.JuMPDict{Uncertain})
    name = printdata(u).name
    for (ind,unc) in zip(keys(u),values(u))
        #if mode == IJuliaMode
        #    uncNames[unc.id] = string(name, "_{", join([string(i) for i in ind],","), "}")
        #else
            uncNames[unc.id] = string(name,  "[", join([string(i) for i in ind],","), "]")
        #end
    end
end
function fill_unc_names(mode, uncNames, u::Array{Uncertain})
    isempty(u) && return
    sizes = size(u)
    m = first(u).m
    rmext = get_robust(m)
    if !haskey(rmext.uncData, u)
        return
    end
    name = rmext.uncData[u].name
    for (ii,unc) in enumerate(u)
        @assert unc.m === m
        ind = ind2sub(sizes, ii)
        #uncNames[unc.id] = if mode === IJuliaMode
        #    string(name, "_{", join(ind, ","), "}")
        #else
        #    string(name,  "[", join(ind, ","), "]")
        #end
        uncNames[unc.id] = string(name,  "[", join(ind, ","), "]")
    end
    return
end

# Handlers to use correct symbols
unc_str(::Type{REPLMode}, u::Uncertain) =
    unc_str(REPLMode, u.m, u.id)
#unc_str(::Type{IJuliaMode}, v::Variable; mathmode=true) =
#    unc_str(IJuliaMode, u.m, u.id, mathmode=mathmode)

#unc_str(::Type{REPLMode}, m::Model, unc::Int) =
#    unc_str(REPLMode, m, unc)
#unc_str(::Type{IJuliaMode}, m::Model, unc::Int; mathmode=true) =
#    math(unc_str(IJuliaMode, m, unc, ijulia_ind_open, ijulia_ind_close), mathmode)


#------------------------------------------------------------------------
## JuMPContainer{Uncertain}
#------------------------------------------------------------------------
Base.print(io::IO, j::Union{JuMPContainer{Uncertain},Array{Uncertain}}) = print(io, cont_str(REPLMode,j))
Base.show( io::IO, j::Union{JuMPContainer{Uncertain},Array{Uncertain}}) = print(io, cont_str(REPLMode,j))
#Base.writemime(io::IO, ::MIME"text/latex", j::JuMPContainer{Uncertain}) =
#    print(io, cont_str(IJuliaMode,j,mathmode=false))
# Generic string converter, called by mode-specific handlers
_getmodel(j::Array{Uncertain}) = first(j).m
_getmodel(j::JuMPContainer) = getmeta(j, :model)
function cont_str(mode, j::Union{JuMPContainer{Uncertain},Array{Uncertain}},
                    sym::PrintSymbols)
    # Check if anything in the container
    if isempty(j)
        name = isa(j, JuMPContainer) ? printdata(j).name : "Empty Array{Uncertain}"
        readline()
        return "$name (no indices)"
    end

    m = _getmodel(j)
    rmext = get_robust(m)
    data = printdata(j)

    # 1. construct the part with uncertain name and indexing
    locvars = map(data.indexexprs) do tmp
        var = tmp.idxvar
        if var == nothing
            return ""
        else
            return string(var)
        end
    end
    num_dims = length(data.indexsets)
    idxvars = Array{String}(num_dims)
    dimidx = 1
    for i in 1:num_dims
        if data.indexexprs[i].idxvar == nothing
            while DIMS[dimidx] in locvars
                dimidx += 1
            end
            if dimidx > length(DIMS)
                error("Unexpectedly ran out of indices")
            end
            idxvars[i] = DIMS[dimidx]
            dimidx += 1
        else
            idxvars[i] = locvars[i]
        end
    end
    name_idx = string(data.name, sym[:ind_open], join(idxvars,","), sym[:ind_close])
    # 2. construct part with what we index over
    idx_sets = sym[:for_all]*" "*join(map(dim->string(idxvars[dim], " ", sym[:in],
                                " ", sym[:open_set],
                                JuMP.cont_str_set(data.indexsets[dim], sym[:dots]),
                                sym[:close_set]), 1:num_dims), ", ")
    # 3. Handle any conditionals
    #if isa(dict, JuMP.JuMPDict) && !isempty(dict.condition)
    #    tail_str *= " s.t. $(join(parse_conditions(j.condition[1]), " and "))"
    #end

    # 4. Bounds and category, if possible, and return final string
    a_var = first(JuMP._values(j))
    unc_cat = rmext.unc_cat[a_var.id]
    unc_lb  = rmext.unc_lower[a_var.id]
    unc_ub  = rmext.unc_upper[a_var.id]
    # Variables may have different bounds, so we can't really print nicely
    # at this time (possibly ever, as they could have been changed post
    # creation, which we'd never be able to handle.
    all_same_lb = true
    all_same_ub = true
    for unc in JuMP._values(j)
        all_same_lb &= rmext.unc_lower[unc.id] == unc_lb
        all_same_ub &= rmext.unc_upper[unc.id] == unc_ub
    end
    str_lb = unc_lb == -Inf ? "-"*sym[:infty] : str_round(unc_lb)
    str_ub = unc_ub == +Inf ?     sym[:infty] : str_round(unc_ub)
    # Special case bounds printing based on the category
    if unc_cat == :Bin  # x in {0,1}
        return "$name_idx $(sym[:in]) $(sym[:open_set])0,1$(sym[:close_set]) $idx_sets"
    end
    # Continuous and Integer
    idx_sets = unc_cat == :Int ? ", $(sym[:integer]), $idx_sets" : " $idx_sets"
    if all_same_lb && all_same_ub
        # Free variable
        unc_lb == -Inf && unc_ub == +Inf && return "$name_idx free$idx_sets"
        # No lower bound
        unc_lb == -Inf && return "$name_idx $(sym[:leq]) $str_ub$idx_sets"
        # No upper bound
        unc_ub == +Inf && return "$name_idx $(sym[:geq]) $str_lb$idx_sets"
        # Range
        return "$str_lb $(sym[:leq]) $name_idx $(sym[:leq]) $str_ub$idx_sets"
    end
    if all_same_lb && !all_same_ub
        unc_lb == -Inf && return "$name_idx $(sym[:leq]) $(sym[:dots])$idx_sets"
        return "$str_lb $(sym[:leq]) $name_idx $(sym[:leq]) $(sym[:dots])$idx_sets"
    end
    if !all_same_lb && all_same_ub
        unc_ub == +Inf && return "$name_idx $(sym[:geq]) $(sym[:dots])$idx_sets"
        return "$(sym[:dots]) $(sym[:leq]) $name_idx $(sym[:leq]) $str_ub$idx_sets"
    end
    return "$(sym[:dots]) $(sym[:leq]) $name_idx $(sym[:leq]) $(sym[:dots])$idx_sets"
end

# Handlers to use correct symbols
cont_str(::Type{REPLMode}, j::JuMPContainer{Uncertain}; mathmode=false) =
    cont_str(REPLMode, j, repl)
#=cont_str(::Type{IJuliaMode}, j::JuMPContainer{Uncertain}; mathmode=true) =
    math(cont_str(IJuliaMode, j, ijulia), mathmode)=#


#------------------------------------------------------------------------
## UncExpr
#------------------------------------------------------------------------
Base.print(io::IO, a::UncExpr) = print(io, aff_str(REPLMode,a))
Base.show( io::IO, a::UncExpr) = print(io, aff_str(REPLMode,a))
#Base.writemime(io::IO, ::MIME"text/latex", a::UncExpr) =
#    print(io, math(aff_str(IJuliaMode,a),false))
# Generic string converter, called by mode-specific handlers
function aff_str(mode, a::UncExpr, show_constant=true)
    # If the expression is empty, return the constant (or 0)
    if length(a.vars) == 0
        return show_constant ? str_round(a.constant) : "0"
    end

    # Get reference to robust part of model
    m  = a.vars[1].m
    rmext = get_robust(m)

    # Collect like terms
    indvec = JuMP.IndexedVector(Float64, rmext.num_uncs)
    for ind in 1:length(a.vars)
        JuMP.addelt!(indvec, a.vars[ind].id, a.coeffs[ind])
    end

    elm = 1
    term_str = Array{String}(2 * length(a.vars))
    # For each non-zero for this model
    for i in 1:indvec.nnz
        idx = indvec.nzidx[i]
        elt = indvec.elts[idx]
        abs(elt) < PRINT_ZERO_TOL && continue  # e.g. x - x

        pre = abs(abs(elt)-1) < PRINT_ZERO_TOL ? "" : str_round(abs(elt)) * " "
        var = unc_str(mode,m,idx)

        term_str[2*elm-1] = elt < 0 ? " - " : " + "
        term_str[2*elm  ] = "$pre$var"
        elm += 1
    end

    if elm == 1
        # Will happen with cancellation of all terms
        # We should just return the constant, if its desired
        return show_constant ? str_round(a.constant) : "0"
    else
        # Correction for very first term - don't want a " + "/" - "
        term_str[1] = (term_str[1] == " - ") ? "-" : ""
        ret = join(term_str[1:2*(elm-1)])
        if abs(a.constant) >= PRINT_ZERO_TOL && show_constant
            ret = string(ret, a.constant < 0 ? " - " : " + ", str_round(abs(a.constant)))
        end
        return ret
    end
end


#------------------------------------------------------------------------
## UncVarExpr
#------------------------------------------------------------------------
Base.print(io::IO, a::UncVarExpr) = print(io, aff_str(REPLMode,a))
Base.show( io::IO, a::UncVarExpr) = print(io, aff_str(REPLMode,a))
#Base.writemime(io::IO, ::MIME"text/latex", a::UncVarExpr) =
#    print(io, math(aff_str(IJuliaMode,a),false))
# Generic string converter, called by mode-specific handlers
# con_str says showConstant = false, because if the constant term is
# just a number for AffExpr. However in our case it might also contain
# uncertains - which we ALWAYS want to show.
# So we "partially" respect it.
function aff_str(mode, a::UncVarExpr, show_constant=true)
    # If no variables, hand off to the constant part
    if length(a.vars) == 0
        return aff_str(mode, a.constant)
    end

    # Get reference to robust part of model
    m  = a.vars[1].m
    rmext = get_robust(m)

    # Don't collect like terms
    term_str = Array{String}(length(a.vars))
    numTerms = 0
    first = true
    for i in 1:length(a.vars)
        numTerms += 1
        uaff = a.coeffs[i]
        varn = getname(a.vars[i])
        prefix = first ? "" : " + "
        # Coefficient expression is a constant
        if length(uaff.vars) == 0
            if abs(uaff.constant) <= PRINT_ZERO_TOL
                # Constant 0 - do not display this term at all
                term_str[numTerms] = ""
            elseif abs(uaff.constant - 1) <= PRINT_ZERO_TOL
                # Constant +1
                term_str[numTerms] = first ? varn : " + $varn"
            elseif abs(uaff.constant + 1) <= PRINT_ZERO_TOL
                # Constant -1
                term_str[numTerms] = first ? "-$varn" : " - $varn"
            else
                # Constant is other than 0, +1, -1
                sign = first ? (uaff.constant < 0 ? "-" : "") :
                               (uaff.constant < 0 ? " - " : " + ")
                term_str[numTerms] = "$(sign)$(str_round(abs(uaff.constant))) $varn"
            end
        # Coefficient expression is a single uncertainty
        elseif length(uaff.vars) == 1
            if abs(uaff.constant) <= PRINT_ZERO_TOL && abs(abs(uaff.coeffs[1]) - 1) <= PRINT_ZERO_TOL
                # No constant, so no (...) needed
                term_str[numTerms] = string(prefix,aff_str(mode,uaff)," ",varn)
            else
                # Constant - need (...)
                ustr = aff_str(mode,uaff)
                # Check its not an all-zero expression
                if ustr == "0"
                    term_str[numTerms] = ""
                else
                    term_str[numTerms] = string(prefix,"(",aff_str(mode,uaff),") ",varn)
                end
            end
        # Coefficient is a more complicated expression
        else
            ustr = aff_str(mode,uaff)
            # Check its not an all-zero expression
            if ustr == "0"
                term_str[numTerms] = ""
            else
                term_str[numTerms] = string(prefix,"(",aff_str(mode,uaff),") ",varn)
            end
        end
        first = false
    end

    # And then connect them up
    ret = join(term_str[1:numTerms], "")

    # Now the constant term
    con_aff = aff_str(REPLMode,a.constant,show_constant)
    if con_aff != "" && con_aff != "0"
        ret = string(ret, " + ", con_aff)
    end

    return ret
end


#------------------------------------------------------------------------
## UncSetNormConstraint
#------------------------------------------------------------------------
Base.print(io::IO, unc::UncSetNormConstraint) = print(io, con_str(REPLMode,unc))
Base.show( io::IO, unc::UncSetNormConstraint) = print(io, con_str(REPLMode,unc))
#Base.writemime(io::IO, ::MIME"text/latex", e::UncSetNormConstraint) =
#    print(io, math(con_str(IJuliaMode,e),false))
# Generic string converter, called by mode-specific handlers
function con_str{P}(mode, unc::UncSetNormConstraint{P})
    normexpr = unc.normexpr
    nrm = normexpr.norm
    cof = normexpr.coeff
    aff = normexpr.aff
    # Coefficient out front
    ret = (cof == 1.0) ? "" : str_round(cof)
    # Norm part
    ret *= "‖" * join(map(t->aff_str(mode,t),nrm.terms),",") * "‖"
    P ==   1 && (ret *= "₁")
    P ==   2 && (ret *= "₂")
    P == Inf && (ret *= "∞")
    # RHS
    ret *= " $(repl[:leq]) "
    @assert length(aff.vars) == 0
    ret *= str_round(-aff.constant)
    return ret
end


#------------------------------------------------------------------------
## Adaptive
#------------------------------------------------------------------------
function adp_str(mode, m::Model, id::Int)
    rd = get_robust(m)
    return rd.adp_names[id] == "" ? "adp_$id" : rd.adp_names[id]
end
adp_str(m::Model, id::Int) = adp_str(REPLMode, m::Model, id::Int)


#------------------------------------------------------------------------
## AdaptExpr
#------------------------------------------------------------------------
Base.show(io::IO, a::AdaptExpr) = print(io, aff_str(REPLMode,a))
# Generic string converter, called by mode-specific handlers
function aff_str(mode, a::AdaptExpr, show_constant=true)
    # If the expression is empty, return the constant (or 0)
    if length(a.vars) == 0
        return show_constant ? str_round(a.constant) : "0"
    end

    # Get reference to model and robust part of model
    m  = a.vars[1].m
    rd = get_robust(m)

    # Collect like terms
    indvec_var = JuMP.IndexedVector(Float64,  m.numCols)
    indvec_adp = JuMP.IndexedVector(Float64, rd.num_adps)
    for ind in 1:length(a.vars)
        if isa(a.vars[ind], Adaptive)
            JuMP.addelt!(indvec_adp, a.vars[ind].id, a.coeffs[ind])
        elseif isa(a.vars[ind], Variable)
            JuMP.addelt!(indvec_var, a.vars[ind].col, a.coeffs[ind])
        end
    end

    elm = 1
    term_str = Array{String}(2 * length(a.vars))
    # For each non-zero
    for i in 1:indvec_var.nnz
        idx = indvec_var.nzidx[i]
        elt = indvec_var.elts[idx]
        abs(elt) < PRINT_ZERO_TOL && continue  # e.g. x - x

        pre = abs(abs(elt)-1) < PRINT_ZERO_TOL ? "" : str_round(abs(elt)) * " "
        var = JuMP.var_str(mode,m,idx)

        term_str[2*elm-1] = elt < 0 ? " - " : " + "
        term_str[2*elm  ] = "$pre$var"
        elm += 1
    end
    for i in 1:indvec_adp.nnz
        idx = indvec_adp.nzidx[i]
        elt = indvec_adp.elts[idx]
        abs(elt) < PRINT_ZERO_TOL && continue  # e.g. x - x

        pre = abs(abs(elt)-1) < PRINT_ZERO_TOL ? "" : str_round(abs(elt)) * " "
        var = adp_str(mode,m,idx)

        term_str[2*elm-1] = elt < 0 ? " - " : " + "
        term_str[2*elm  ] = "$pre$var"
        elm += 1
    end

    if elm == 1
        # Will happen with cancellation of all terms
        # We should just return the constant, if its desired
        return show_constant ? str_round(a.constant) : "0"
    else
        # Correction for very first term - don't want a " + "/" - "
        term_str[1] = (term_str[1] == " - ") ? "-" : ""
        ret = join(term_str[1:2*(elm-1)])
        if abs(a.constant) >= PRINT_ZERO_TOL && show_constant
            ret = string(ret, a.constant < 0 ? " - " : " + ", str_round(abs(a.constant)))
        end
        return ret
    end
end
