#############################################################################
# JuMPeR
# Julia for Mathematical Programming - extension for Robust Optimization
# See http://github.com/IainNZ/JuMPeR.jl
#############################################################################

#############################################################################
#### Robust model-specific printing
#############################################################################

# checkNameStatus
# [Not exported] Initially uncertains defined as indexed with @defUnc do not
# have names because generating them is relatively slow compared to just
# creating the uncertain. We lazily generate them only when someone requires
# them. This function checks if we have done that already, and if not,
# generates them
function checkUncNameStatus(rm::Model)
    rd = getRobust(rm)
    for i in 1:rd.numUncs
        if rd.uncNames[i] == ""
            fillUncNames(rm)
            return
        end
    end
end

# fillUncNames
# Go through every unc JuMPDict and create names for the uncertainties
# that go with every dictionary.
function fillUncNames(rm::Model)
    rd = getRobust(rm)
    for dict in rd.dictList
        idxsets = dict.indexsets
        lengths = map(length, idxsets)
        N = length(idxsets)
        name = dict.name
        cprod = cumprod([lengths...])
        for (ind,unc) in enumerate(dict.innerArray)
            setName(unc,string("$name[$(idxsets[1][mod1(ind,lengths[1])])", [ ",$(idxsets[i][int(ceil(mod1(ind,cprod[i]) / cprod[i-1]))])" for i=2:N ]..., "]"))
        end
    end
end


function printRobust(m::Model)
    checkUncNameStatus(m)
    rd = getRobust(m)
    # First, display normal model stuff
    print(m)
    println("Uncertain constraints:")
    for c in rd.uncertainconstr
        println(conToStr(c))
    end
    println("Uncertainty set:")
    for uc in rd.uncertaintyset
        println(conToStr(uc))
    end
    for elc in rd.normconstraints
        printEll(m,elc)
    end
    for unc in 1:rd.numUncs
        println("$(rd.uncLower[unc]) <= $(rd.uncNames[unc]) <= $(rd.uncUpper[unc])")
    end
end


#############################################################################
#############################################################################

function affToStr(a::UAffExpr, showConstant=true)
    # If the expression is empty, return the constant (or 0)
    if length(a.vars) == 0
        return showConstant ? str_round(a.constant) : "0"
    end

    # Get reference to robust part of model
    robdata = getRobust(a.vars[1].m)
    checkUncNameStatus(a.vars[1].m)

    # Collect like terms
    indvec = IndexedVector(Float64, robdata.numUncs)
    for ind in 1:length(a.vars)
        addelt!(indvec, a.vars[ind].unc, a.coeffs[ind])
    end

    # Stringify the terms
    elm = 1
    term_str = Array(UTF8String, 2*length(a.vars))
    for i in 1:indvec.nnz
        idx = indvec.nzidx[i]
        elt = indvec.elts[idx]
        abs(elt) < JuMP.PRINT_ZERO_TOL && continue  # e.g. x - x

        pre = abs(abs(elt)-1) < JuMP.PRINT_ZERO_TOL ? "" : str_round(abs(elt)) * " "
        var = robdata.uncNames[idx]

        term_str[2*elm-1] = elt < 0 ? " - " : " + "
        term_str[2*elm  ] = "$pre$var"
        elm += 1
    end
        
    if elm == 1
        # Will happen with cancellation of all terms
        # We should just return the constant, if its desired
        return showConstant ? str_round(a.constant) : "0"
    else
        # Correction for very first term - don't want a " + "/" - "
        term_str[1] = (term_str[1] == " - ") ? "-" : ""
        ret = join(term_str[1:2*(elm-1)])
        if abs(a.constant) >= JuMP.PRINT_ZERO_TOL && showConstant
            ret = string(ret, a.constant < 0 ? " - " : " + ", str_round(abs(a.constant)))
        end
        return ret
    end
end

#############################################################################
#############################################################################

# Now conToStr says showConstant = false, because if the constant term is
# just a number for AffExpr. However in our case it might also contain
# uncertains - which we ALWAYS want to show. So we "partially" respect it.
function affToStr(a::FullAffExpr, showConstant=true)
    const ZEROTOL = 1e-20

    # If no variables, hand off to the constant part
    if length(a.vars) == 0
        # Will do the right thing even with showContant=true or false
        return affToStr(a.constant)
    end

    # Get reference to robust part of model
    robdata = getRobust(a.vars[1].m)
    checkUncNameStatus(a.vars[1].m)

    # Stringify the terms - we don't collect like terms
    termStrings = Array(UTF8String, length(a.vars))
    numTerms = 0
    first = true
    for i in 1:length(a.vars)
        numTerms += 1
        uaff = a.coeffs[i]
        varn = getName(a.vars[i])
        prefix = first ? "" : " + "
        # Coefficient expression is a constant
        if length(uaff.vars) == 0
            if abs(uaff.constant) <= ZEROTOL
                # Constant 0 - do not display this term at all
                termStrings[numTerms] = ""
            elseif abs(uaff.constant - 1) <= ZEROTOL
                # Constant +1
                termStrings[numTerms] = first ? varn : " + $varn"
            elseif abs(uaff.constant + 1) <= ZEROTOL
                # Constant -1
                termStrings[numTerms] = first ? "-$varn" : " - $varn"
            else
                # Constant is other than 0, +1, -1 
                if first
                    sign = uaff.constant < 0 ? "-" : ""
                    termStrings[numTerms] = "$sign$(str_round(abs(uaff.constant))) $varn"
                else
                    sign = uaff.constant < 0 ? "-" : "+"
                    termStrings[numTerms] = " $sign $(str_round(abs(uaff.constant))) $varn"
                end
            end
        # Coefficient expression is a single uncertainty
        elseif length(uaff.vars) == 1
            if abs(uaff.constant) <= ZEROTOL && abs(abs(uaff.coeffs[1]) - 1) <= ZEROTOL
                # No constant, so no (...) needed
                termStrings[numTerms] = string(prefix,affToStr(uaff)," ",varn)
            else
                # Constant - need (...)
                termStrings[numTerms] = string(prefix,"(",affToStr(uaff),") ",varn)
            end
        # Coefficient is a more complicated expression
        else
            termStrings[numTerms] = string(prefix,"(",affToStr(uaff),") ",varn)
        end
        first = false
    end

    # And then connect them up with +s
    ret = join(termStrings[1:numTerms], "")
    
    # Now the constant term
    con_aff = affToStr(a.constant,showConstant)
    if con_aff != "" && con_aff != "0"
        ret = string(ret, " + ", con_aff)
    end

    return ret
end