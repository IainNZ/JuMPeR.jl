macro defUnc(m, x, extra...)
    m = esc(m)
    if isexpr(x,:comparison)
        # we have some bounds
        if x.args[2] == :>=
            if length(x.args) == 5
                error("Use the form lb <= unc <= ub instead of ub >= unc >= lb")
            end
            @assert length(x.args) == 3
            # lower bounds, no upper
            lb = esc(x.args[3])
            ub = Inf
            var = x.args[1]
        elseif x.args[2] == :<=
            if length(x.args) == 5
                # lb <= x <= u
                lb = esc(x.args[1])
                if (x.args[4] != :<=)
                    error("Expected <= operator")
                end
                ub = esc(x.args[5])
                var = x.args[3]
            else
                # x <= u
                ub = esc(x.args[3])
                lb = -Inf
                var = x.args[1]
            end
        end
    else
        var = x
        lb = -Inf
        ub = Inf
    end
    
    t = :Cont
    if length(extra) > 0
        gottype = 0
        if extra[1] == :Int || extra[1] == :Bin
            gottype = 1
            if extra[1] == :Int
                t = :Int
            else
                if lb != -Inf || ub != Inf
                    error("Bounds may not be specified for binary variables. These are always taken to have a lower bound of 0 and upper bound of 1.")
                end
                t = :Int
                lb = 0.0
                ub = 1.0
            end
        end
    end


    if isa(var,Symbol)
        # easy case
        return quote
            $(esc(var)) = Uncertain($m,$lb,$ub,$(Meta.quot(t)),$(string(var)))
            nothing
        end
    else
        if !isexpr(var,:ref)
            error("Syntax error: Expected $var to be of form var[...]")
        end
        varname = esc(var.args[1])
        idxvars = {}
        idxsets = {}
        refcall = Expr(:ref,varname)
        for s in var.args[2:end]
            if isa(s,Expr) && s.head == :(=)
                idxvar = esc(s.args[1])
                idxset = esc(s.args[2])
            else
                idxvar = gensym()
                idxset = esc(s)
            end
            push!(idxvars, idxvar)
            push!(idxsets, idxset)
            push!(refcall.args, idxvar)
        end
        code = :( $(refcall) = Uncertain($m, $lb, $ub, $(Meta.quot(t))) )
        for (idxvar, idxset) in zip(reverse(idxvars),reverse(idxsets))
            code = quote
                for $idxvar in $idxset
                    $code
                end
            end
        end
        
        mac = Expr(:macrocall,symbol("@gendict"),varname,:Uncertain,idxsets...)
        addDict = :( push!($(m).ext[:Robust].dictList, $varname) )
        code = quote 
            $mac
            $code
            $addDict
            nothing
        end
        return code
    end
end

# Stuff to make JuMP macros work with Uncertains
(*)(u::Uncertain) = u
function addToExpression(aff::GenericAffExpr,c,x)
    return aff + c*x
end
