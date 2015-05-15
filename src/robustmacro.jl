macro defUnc(m, x, extra...)
    #-------- DIRECT COPY BEGIN ------- (except 'variable' -> 'uncertain')
    m = esc(m)
    # Identify the variable bounds. Four (legal) possibilities are "x >= lb",
    # "x <= ub", "lb <= x <= ub", or just plain "x"
    if isexpr(x,:comparison)
        # We have some bounds
        if x.args[2] == :>= || x.args[2] == :≥
            if length(x.args) == 5
                # ub >= x >= lb
                x.args[4] == :>= || x.args[4] == :≥ || error("Invalid uncertain bounds")
                var = x.args[3]
                lb = esc(x.args[5])
                ub = esc(x.args[1])
            else
                # x >= lb
                var = x.args[1]
                @assert length(x.args) == 3
                lb = esc(x.args[3])
                ub = Inf
            end
        elseif x.args[2] == :<= || x.args[2] == :≤
            if length(x.args) == 5
                # lb <= x <= u
                var = x.args[3]
                (x.args[4] != :<= && x.args[4] != :≤) &&
                    error("in @defUnc ($var): expected <= operator after uncertain name.")
                lb = esc(x.args[1])
                ub = esc(x.args[5])
            else
                # x <= ub
                var = x.args[1]
                # NB: May also be lb <= x, which we do not support
                #     We handle this later in the macro
                @assert length(x.args) == 3
                ub = esc(x.args[3])
                lb = -Inf
            end
        else
            # Its a comparsion, but not using <= ... <=
            error("in @defUnc ($(string(x))): use the form lb <= ... <= ub.")
        end
    else
        # No bounds provided - free variable
        # If it isn't, e.g. something odd like f(x), we'll handle later
        var = x
        lb = -Inf
        ub = Inf
    end
    #-------- DIRECT COPY END ------- 
    
    #-------- DIRECT COPY BEGIN ------- 
    # - except 'variable' -> 'uncertain'
    # - except no semicont, semiint
    # - except column format
    t = :Cont
    gottype = 0
    if length(extra) > 0
        if extra[1] in [:Bin, :Int]
            gottype = 1
            t = extra[1]
        end

        if t == :Bin 
            if (lb != -Inf || ub != Inf) && !(lb == 0.0 && ub == 1.0)
            error("in @defUnc ($var): bounds other than [0, 1] may not be specified for binary uncertains.\nThese are always taken to have a lower bound of 0 and upper bound of 1.")
            else
                lb = 0.0
                ub = 1.0
            end
        end

        gottype == 0 &&
            error("in @defUnc ($var): syntax error")
    end
    #-------- DIRECT COPY END ------- 

    #-------- DIRECT COPY BEGIN ------- 
    # - except Variable -> Uncertain
    if isa(var,Symbol)
        # Easy case - a single uncertain
        return quote
            $(esc(var)) = Uncertain($m,$lb,$ub,$(Meta.quot(t)),$(string(var)))
            nothing
        end
    end
    @assert isa(var,Expr)
    #-------- DIRECT COPY END ------- 


    #-------- DIRECT COPY BEGIN ------- 
    # - Variable -> Uncertain
    # - m.dictlist -> m.ext[:Robust].dictList
    refcall, idxvars, idxsets, idxpairs = JuMP.buildrefsets(var)
    code = :( $(refcall) = Uncertain($m, $lb, $ub, $(Meta.quot(t))) )
    looped = JuMP.getloopedcode(var, code, :(), idxvars, idxsets, idxpairs, :Uncertain)
    varname = esc(JuMP.getname(var))
    return quote 
        $looped
        push!($(m).ext[:Robust].dictList, $varname)
        $varname
    end
end

# Stuff to make JuMP macros work with Uncertains - should probably be
# typed tighter, but seems to work OK.
(*)(u::Uncertain) = u
function JuMP.addToExpression(aff::GenericAffExpr, c::Number, x::Union(Number,UAffExpr))
    return aff + c*x
end
function JuMP.addToExpression(aff::GenericAffExpr, c::Union(Number,AffExpr), x::Uncertain)
    return aff + c*x
end
function JuMP.addToExpression(aff::GenericAffExpr, c::Union(Number,UAffExpr), x::Variable)
    return aff + c*x
end
function JuMP.addToExpression(val::Real, c::Number, x::Union(UAffExpr,FullAffExpr))
    return val + c*x
end

function JuMP._construct_constraint!(faff::FullAffExpr, sense::Symbol)
    JuMP._canonicalize_sense(sense)
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