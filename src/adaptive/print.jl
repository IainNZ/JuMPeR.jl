#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2015: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# src/adaptive/print.jl
# Adaptive robust optimization support - printing
#-----------------------------------------------------------------------

function adp_str(mode, m::Model, id::Int)
    rd = getRobust(m)
    return rd.adpNames[id] == "" ? "adp_$id" : rd.adpNames[id]
end
adp_str(m::Model, id::Int) = adp_str(REPLMode, m::Model, id::Int)


#------------------------------------------------------------------------
## AdaptAffExpr
#------------------------------------------------------------------------
Base.print(io::IO, a::AdaptAffExpr) = print(io, aff_str(REPLMode,a))
Base.show( io::IO, a::AdaptAffExpr) = print(io, aff_str(REPLMode,a))
#Base.writemime(io::IO, ::MIME"text/latex", a::AdaptAffExpr) =
#    print(io, math(aff_str(IJuliaMode,a),false))
# Generic string converter, called by mode-specific handlers
function aff_str(mode, a::AdaptAffExpr, show_constant=true)
    # If the expression is empty, return the constant (or 0)
    if length(a.vars) == 0
        return show_constant ? str_round(a.constant) : "0"
    end

    # Get reference to model and robust part of model
    m  = a.vars[1].m
    rd = getRobust(m)

    # Collect like terms
    indvec_var = IndexedVector(Float64,  m.numCols)
    indvec_adp = IndexedVector(Float64, rd.numAdapt)
    for ind in 1:length(a.vars)
        if isa(a.vars[ind], Adaptive)
            addelt!(indvec_adp, a.vars[ind].id, a.coeffs[ind])
        elseif isa(a.vars[ind], Variable)
            addelt!(indvec_var, a.vars[ind].col, a.coeffs[ind])
        end
    end

    elm = 1
    term_str = Array(UTF8String, 2*length(a.vars))
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

# Backwards compatability shim
affToStr(a::AdaptAffExpr) = aff_str(REPLMode,a)