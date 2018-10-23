#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2016: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# src/robustmacro.jl
# This file attempts to use as much of the machinery provided by JuMP
# to provide the @uncertain functionality for Uncertains, just like @variable
# does for variables. @uncertain is essentially a direct copy of @variable,
# with the differences clearly marked. If something is not marked,
# assume it is a direct copy.
#-----------------------------------------------------------------------

import JuMP: assert_validmodel, validmodel, esc_nonconstant
import JuMP: getloopedcode, buildrefsets
import JuMP: storecontainerdata, isdependent, JuMPArray, JuMPContainerData, pushmeta!
using Base.Meta

storecontainerdata_unc(m::Model, variable, varname, idxsets, idxpairs, condition) =
    get_robust(m).uncData[variable] = JuMPContainerData(varname, idxsets, idxpairs, condition)

macro uncertain(args...)
    # Synchronized to JuMP's @variable at commit:
    # 6d888fe625138074c71322e2ebdf9cf14ed3444c (April 27 2016)
    # MODIFICATION: error message
    length(args) <= 1 &&
        error("in @uncertain: expected model as first argument, then uncertain parameter information.")
    m = esc(args[1])
    x = args[2]
    extra = vcat(args[3:end]...)

    t = :Cont
    gottype = false
    haslb = false
    hasub = false
    # Identify the variable bounds. Five (legal) possibilities are "x >= lb",
    # "x <= ub", "lb <= x <= ub", "x == val", or just plain "x"
    if VERSION < v"0.5.0-dev+3231"
        x = JuMP.comparison_to_call(x)
    end
    if isexpr(x,:comparison) # two-sided
        haslb = true
        hasub = true
        if x.args[2] == :>= || x.args[2] == :≥
            # ub >= x >= lb
            # MODIFICATION: error message
            x.args[4] == :>= || x.args[4] == :≥ || error("Invalid uncertain parameter bounds")
            var = x.args[3]
            lb = esc_nonconstant(x.args[5])
            ub = esc_nonconstant(x.args[1])
        elseif x.args[2] == :<= || x.args[2] == :≤
            # lb <= x <= u
            var = x.args[3]
            # MODIFICATION: error message
            (x.args[4] != :<= && x.args[4] != :≤) &&
                error("in @uncertain ($var): expected <= operator after uncertain parameter name.")
            lb = esc_nonconstant(x.args[1])
            ub = esc_nonconstant(x.args[5])
        else
            # MODIFICATION: error message
            error("in @uncertain ($(string(x))): use the form lb <= ... <= ub.")
        end
    elseif isexpr(x,:call)
        if x.args[1] == :>= || x.args[1] == :≥
            # x >= lb
            var = x.args[2]
            @assert length(x.args) == 3
            lb = esc_nonconstant(x.args[3])
            haslb = true
            ub = Inf
        elseif x.args[1] == :<= || x.args[1] == :≤
            # x <= ub
            var = x.args[2]
            # NB: May also be lb <= x, which we do not support
            #     We handle this later in the macro
            @assert length(x.args) == 3
            ub = esc_nonconstant(x.args[3])
            hasub = true
            lb = -Inf
        elseif x.args[1] == :(==)
            # fixed variable
            var = x.args[2]
            @assert length(x.args) == 3
            lb = esc(x.args[3])
            haslb = true
            ub = esc(x.args[3])
            hasub = true
            gottype = true
            # MODIFICATION: no :Fixed type for uncertain parameters
            t = :Cont
        else
            # Its a comparsion, but not using <= ... <=
            # MODIFICATION: error message
            error("in @uncertain: unexpected syntax $(string(x)).")
        end
    else
        # No bounds provided - free variable
        # If it isn't, e.g. something odd like f(x), we'll handle later
        var = x
        lb = -Inf
        ub = Inf
    end

    # separate out keyword arguments
    kwargs = filter(ex->isexpr(ex,:kw), extra)
    extra = filter(ex->!isexpr(ex,:kw), extra)

    # process keyword arguments
    # MODIFICATION: no column generation, no initial values
    quotvarname = quot(getname(var))
    escvarname  = esc(getname(var))
    symvarname = Symbol(getname(var))
    for ex in kwargs
        kwarg = ex.args[1]
        if kwarg == :basename
            quotvarname = esc(ex.args[2])
        elseif kwarg == :lowerbound
            # MODIFICATION: error message
            haslb && error("Cannot specify uncertain parameter lowerbound twice")
            lb = esc_nonconstant(ex.args[2])
            haslb = true
        elseif kwarg == :upperbound
            # MODIFICATION: error message
            hasub && error("Cannot specify uncertain parameter upperbound twice")
            ub = esc_nonconstant(ex.args[2])
            hasub = true
        else
            # MODIFICATION: error message
            error("in @uncertain ($var): Unrecognized keyword argument $kwarg")
        end
    end

    sdp = any(t -> (t == :SDP), extra)
    symmetric = (sdp || any(t -> (t == :Symmetric), extra))
    extra = filter(x -> (x != :SDP && x != :Symmetric), extra) # filter out SDP and sym tag
    # MODIFICATION: no SDP
    if sdp || symmetric
        error("in @uncertain ($var): SDP and Symmetric not supported.")
    end

    # Determine variable type (if present).
    # Types: default is continuous (reals)
    if length(extra) > 0
        # MODIFICATION: no fixed
        if extra[1] in [:Bin, :Int, :SemiCont, :SemiInt]
            gottype = true
            t = extra[1]
        end

        if t == :Bin
            if (lb != -Inf || ub != Inf) && !(lb == 0.0 && ub == 1.0)
            error("in @uncertain ($var): bounds other than [0, 1] may not be specified for binary uncertain parameters.\nThese are always taken to have a lower bound of 0 and upper bound of 1.")
            else
                lb = 0.0
                ub = 1.0
            end
        end
        # MODIFICATION: error message
        !gottype && error("in @variable ($var): syntax error")
    end

    # Handle the column generation functionality
    # MODIFICATION: skipped

    if isa(var,Symbol)
        # Easy case - a single variable
        # MODIFICATION: skip SDP check
        return assert_validmodel(m, quote
            # MODIFICATION: Variable to Uncertain, no start value
            $(esc(var)) = Uncertain($m,$lb,$ub,$(quot(t)),string($quotvarname))
            # MODIFICATION: No equivalent to registering variables
        end)
    end
    # MODIFICATION: error message
    isa(var,Expr) || error("in @uncertain: expected $var to be an uncertain parameter name")

    # We now build the code to generate the variables (and possibly the JuMPDict
    # to contain them)
    refcall, idxvars, idxsets, idxpairs, condition = buildrefsets(var)
    clear_dependencies(i) = (isdependent(idxvars,idxsets[i],i) ? nothing : idxsets[i])

    # MODIFICATION: Variable to Uncertain, no start value
    code = :( $(refcall) = Uncertain($m, $lb, $ub, $(quot(t)), JuMP.EMPTYSTRING) )
    # MODIFICATION: No symmetric support
    # if symmetric
    #     ...
    # else

    # MODIFICATION: Variable to Uncertain
    looped = getloopedcode(getname(var), code, condition, idxvars, idxsets, idxpairs, :Uncertain)
    return assert_validmodel(m, quote
        $looped
        # MODIFICATION: use RMExt dictList instead of model dictList
        push!(get_robust($(m)).dictList, $symvarname)
        # MODIFICATION: No equivalent to registering variables
        # MODIFICATION: storecontainerdata -> storecontainerdata_unc
        storecontainerdata_unc($m, $symvarname, $quotvarname,
                           $(Expr(:tuple,map(clear_dependencies,1:length(idxsets))...)),
                           $idxpairs, $(quot(condition)))
        isa($symvarname, JuMPContainer) && pushmeta!($symvarname, :model, $m)
        $escvarname = $symvarname
    end)
end


macro adaptive(args...)
    # Synchronized to JuMP's @variable at commit:
    # 6d888fe625138074c71322e2ebdf9cf14ed3444c (April 27 2016)
    # MODIFICATION: error message
    length(args) <= 1 &&
        error("in @adaptive: expected model as first argument, then variable information.")
    m = esc(args[1])
    x = args[2]
    extra = vcat(args[3:end]...)

    t = :Cont
    gottype = false
    haslb = false
    hasub = false
    # Identify the variable bounds. Five (legal) possibilities are "x >= lb",
    # "x <= ub", "lb <= x <= ub", "x == val", or just plain "x"
    if VERSION < v"0.5.0-dev+3231"
        x = JuMP.comparison_to_call(x)
    end
    if isexpr(x,:comparison) # two-sided
        haslb = true
        hasub = true
        if x.args[2] == :>= || x.args[2] == :≥
            # ub >= x >= lb
            # MODIFICATION: error message
            x.args[4] == :>= || x.args[4] == :≥ || error("Invalid adaptive variable bounds")
            var = x.args[3]
            lb = esc_nonconstant(x.args[5])
            ub = esc_nonconstant(x.args[1])
        elseif x.args[2] == :<= || x.args[2] == :≤
            # lb <= x <= u
            var = x.args[3]
            # MODIFICATION: error message
            (x.args[4] != :<= && x.args[4] != :≤) &&
                error("in @adaptive ($var): expected <= operator after adaptive variable name.")
            lb = esc_nonconstant(x.args[1])
            ub = esc_nonconstant(x.args[5])
        else
            # MODIFICATION: error message
            error("in @adaptive ($(string(x))): use the form lb <= ... <= ub.")
        end
    elseif isexpr(x,:call)
        if x.args[1] == :>= || x.args[1] == :≥
            # x >= lb
            var = x.args[2]
            @assert length(x.args) == 3
            lb = esc_nonconstant(x.args[3])
            haslb = true
            ub = Inf
        elseif x.args[1] == :<= || x.args[1] == :≤
            # x <= ub
            var = x.args[2]
            # NB: May also be lb <= x, which we do not support
            #     We handle this later in the macro
            @assert length(x.args) == 3
            ub = esc_nonconstant(x.args[3])
            hasub = true
            lb = -Inf
        elseif x.args[1] == :(==)
            # fixed variable
            var = x.args[2]
            @assert length(x.args) == 3
            lb = esc(x.args[3])
            haslb = true
            ub = esc(x.args[3])
            hasub = true
            gottype = true
            # MODIFICATION: no :Fixed type for uncertain parameters
            t = :Cont
        else
            # Its a comparsion, but not using <= ... <=
            # MODIFICATION: error message
            error("in @adaptive: unexpected syntax $(string(x)).")
        end
    else
        # No bounds provided - free variable
        # If it isn't, e.g. something odd like f(x), we'll handle later
        var = x
        lb = -Inf
        ub = Inf
    end

    # separate out keyword arguments
    println("================")
    show_sexpr(extra)
    println("================")
    kwargs = filter(ex->isexpr(ex,:kw), extra)
    extra = filter(ex->!isexpr(ex,:kw), extra)

    # process keyword arguments
    # MODIFICATION: no column generation, no initial values
    quotvarname = quot(getname(var))
    escvarname  = esc(getname(var))
    symvarname = Symbol(getname(var))
    # MODIFICATION: adaptive-only things
    policy = :Static
    depends_on = Uncertain[]
    for ex in kwargs
        kwarg = ex.args[1]
        if kwarg == :basename
            quotvarname = esc(ex.args[2])
        elseif kwarg == :lowerbound
            # MODIFICATION: error message
            haslb && error("Cannot specify adaptive variable lowerbound twice")
            lb = esc_nonconstant(ex.args[2])
            haslb = true
        elseif kwarg == :upperbound
            # MODIFICATION: error message
            hasub && error("Cannot specify adaptive variable upperbound twice")
            ub = esc_nonconstant(ex.args[2])
            hasub = true
        elseif kwarg == :policy
            policy = ex.args[2]
        elseif kwarg == :depends_on
            depends_on = esc(ex.args[2])
        else
            # MODIFICATION: error message
            error("in @adaptive ($var): Unrecognized keyword argument $kwarg")
        end
    end

    sdp = any(t -> (t == :SDP), extra)
    symmetric = (sdp || any(t -> (t == :Symmetric), extra))
    extra = filter(x -> (x != :SDP && x != :Symmetric), extra) # filter out SDP and sym tag
    # MODIFICATION: no SDP
    if sdp || symmetric
        error("in @adaptive ($var): SDP and Symmetric not supported.")
    end

    # Determine variable type (if present).
    # Types: default is continuous (reals)
    if length(extra) > 0
        # MODIFICATION: no fixed
        if extra[1] in [:Bin, :Int, :SemiCont, :SemiInt]
            gottype = true
            t = extra[1]
        end

        if t == :Bin
            if (lb != -Inf || ub != Inf) && !(lb == 0.0 && ub == 1.0)
            error("in @adaptive ($var): bounds other than [0, 1] may not be specified for binary adaptive variables.\nThese are always taken to have a lower bound of 0 and upper bound of 1.")
            else
                lb = 0.0
                ub = 1.0
            end
        end
        # MODIFICATION: error message
        !gottype && error("in @adaptive ($var): syntax error")
    end

    # Handle the column generation functionality
    # MODIFICATION: skipped

    if isa(var,Symbol)
        # Easy case - a single variable
        # MODIFICATION: skip SDP check
        return assert_validmodel(m, quote
            # MODIFICATION: Variable to Adaptive, no start value, etc
            $(esc(var)) = Adaptive($m,$lb,$ub,$(quot(t)),string($quotvarname),
                                    $(quot(policy)), $(depends_on))
            # MODIFICATION: No equivalent to registering variables
        end)
    end
    # MODIFICATION: error message
    isa(var,Expr) || error("in @adaptive: expected $var to be an adaptive variable name")

    # We now build the code to generate the variables (and possibly the JuMPDict
    # to contain them)
    refcall, idxvars, idxsets, idxpairs, condition = buildrefsets(var)
    clear_dependencies(i) = (isdependent(idxvars,idxsets[i],i) ? nothing : idxsets[i])

    # MODIFICATION: Variable to Adaptive, no start value
    code = :( $(refcall) = Adaptive($m, $lb, $ub, $(quot(t)), JuMP.EMPTYSTRING,
                                    $(quot(policy)), $(depends_on)) )
    # MODIFICATION: No symmetric support
    # if symmetric
    #     ...
    # else
    # MODIFICATION: Variable to Adaptive
    looped = getloopedcode(getname(var), code, condition, idxvars, idxsets, idxpairs, :Adaptive)
    return assert_validmodel(m, quote
        $looped
        # MODIFICATION: no fancy name stuff
        $escvarname = $symvarname
    end)
end


function JuMP.constructconstraint!(faff::UncVarExpr, sense::Symbol)
    offset = faff.constant.constant
    faff.constant.constant = 0.0
    if sense == :(<=) || sense == :≤
        return UncConstraint(faff, -Inf, -offset)
    elseif sense == :(>=) || sense == :≥
        return UncConstraint(faff, -offset, Inf)
    elseif sense == :(==)
        return UncConstraint(faff, -offset, -offset)
    else
        error("Cannot handle ranged constraint")
    end
end


function JuMP.constructconstraint!{P}(
    normexpr::GenericNormExpr{P,Float64,Uncertain}, sense::Symbol)
    if sense == :(<=)
        UncSetNormConstraint( normexpr)
    elseif sense == :(>=)
        UncSetNormConstraint(-normexpr)
    else
        error("Invalid sense $sense in norm constraint")
    end
end


function JuMP.constructconstraint!(vaff::AdaptExpr, sense::Symbol)
    offset = vaff.constant
    vaff.constant = 0.0
    if sense == :(<=) || sense == :≤
        return AdaptConstraint(vaff, -Inf, -offset)
    elseif sense == :(>=) || sense == :≥
        return AdaptConstraint(vaff, -offset, Inf)
    elseif sense == :(==)
        return AdaptConstraint(vaff, -offset, -offset)
    else
        error("Cannot handle ranged constraint")
    end
end
