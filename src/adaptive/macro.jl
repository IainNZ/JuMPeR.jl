#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2015: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# src/adaptive/macro.jl
# Adaptive robust optimization support - macros
#-----------------------------------------------------------------------

export @defAdaptVar

macro defAdaptVar(args...)
    length(args) <= 1 &&
        error("in @defAdaptVar: expected model as first argument, then variable information.")
    m = esc(args[1])
    x = args[2]
    extra = vcat(args[3:end]...)

    # Identify the variable bounds. Five (legal) possibilities are "x >= lb",
    # "x <= ub", "lb <= x <= ub", or just plain "x"
    if isexpr(x,:comparison)
        # We have some bounds
        if x.args[2] == :>= || x.args[2] == :≥
            if length(x.args) == 5
                # ub >= x >= lb
                x.args[4] == :>= || x.args[4] == :≥ || error("Invalid variable bounds")
                var = x.args[3]
                lb = esc_nonconstant(x.args[5])
                ub = esc_nonconstant(x.args[1])
            else
                # x >= lb
                var = x.args[1]
                @assert length(x.args) == 3
                lb = esc_nonconstant(x.args[3])
                ub = Inf
            end
        elseif x.args[2] == :<= || x.args[2] == :≤
            if length(x.args) == 5
                # lb <= x <= u
                var = x.args[3]
                (x.args[4] != :<= && x.args[4] != :≤) &&
                    error("in @defAdaptVar ($var): expected <= operator after variable name.")
                lb = esc_nonconstant(x.args[1])
                ub = esc_nonconstant(x.args[5])
            else
                # x <= ub
                var = x.args[1]
                # NB: May also be lb <= x, which we do not support
                #     We handle this later in the macro
                @assert length(x.args) == 3
                ub = esc_nonconstant(x.args[3])
                lb = -Inf
            end
        else
            # Its a comparsion, but not using <= ... <=
            error("in @defAdaptVar ($(string(x))): use the form lb <= ... <= ub.")
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
    varcat = :Cont
    policy = :Static
    stage  = 0 
    depends_on = Uncertain[]
    for ex in kwargs
        if ex.args[1] == :policy
            policy = ex.args[2]
        elseif ex.args[1] == :stage
            stage = esc(ex.args[2])
        elseif ex.args[1] == :depends_on
            depends_on = esc(ex.args[2])
        else
            error("in @defAdaptVar ($var): Unrecognized keyword argument $(ex.args[1])")
        end
    end

    # Determine variable type (if present).
    # Types: default is continuous (reals)
    if length(extra) > 0
        if extra[1] in [:Bin, :Int]
            gottype = true
            varcat = extra[1]
        end

        if t == :Bin
            if (lb != -Inf || ub != Inf) && !(lb == 0.0 && ub == 1.0)
            error("in @defAdaptVar ($var): bounds other than [0, 1] may not be specified for binary variables.\nThese are always taken to have a lower bound of 0 and upper bound of 1.")
            else
                lb = 0.0
                ub = 1.0
            end
        end

        !gottype && error("in @defAdaptVar ($var): syntax error")
    end

    if isa(var,Symbol)
        # Easy case - a single variable
        return assert_validmodel(m, quote
            $(esc(var)) = 
                Adaptive($m, $(utf8(string(var))),
                                    $lb, $ub, $(quot(varcat)),
                                    $(quot(policy)), $stage, $(depends_on))
            #registervar($m, $(quot(var)), $(esc(var)))
        end)
    end
    isa(var,Expr) || error("in @defAdaptVar: expected $var to be a variable name")

    
    # We now build the code to generate the variables (and possibly the JuMPDict
    # to contain them)
    refcall, idxvars, idxsets, idxpairs, condition = buildrefsets(var)
    clear_dependencies(i) = (isdependent(idxvars,idxsets[i],i) ? nothing : idxsets[i])
    code = :( $(refcall) = Adaptive($m, "", $lb, $ub,
                                        $(quot(varcat)), $(quot(policy)),
                                        $stage, $(depends_on)) )
    looped = getloopedcode(var, code, condition, idxvars, idxsets, idxpairs, :Adaptive)
    varname = esc(getname(var))
    return assert_validmodel(m, quote
        $looped
        #push!($(m).dictList, $varname)
        #registervar($m, $(quot(getname(var))), $varname)
        #storecontainerdata($m, $varname, $(quot(getname(var))),
        #                   $(Expr(:tuple,map(clear_dependencies,1:length(idxsets))...)),
        #                   $idxpairs, $(quot(condition)))
        #isa($varname, JuMPContainer) && pushmeta!($varname, :model, $m)
        $varname
    end)

end


function JuMP.constructconstraint!(vaff::AdaptAffExpr, sense::Symbol)
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