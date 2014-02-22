macro defUnc(m, x, extra...)
    m = esc(m)
    if isexpr(x,:comparison)
        # we have some bounds
        if x.args[2] == :>=
            if length(x.args) == 5
                error("Use the form lb <= var <= ub instead of ub >= var >= lb")
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
    
    if isa(var,Symbol)
        # easy case
        return quote
            $(esc(var)) = Uncertain($m,$lb,$ub,$(string(var)))
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
        code = :( $(refcall) = Uncertain($m, $lb, $ub) )
        for (idxvar, idxset) in zip(reverse(idxvars),reverse(idxsets))
            code = quote
                for $idxvar in $idxset
                    $code
                end
            end
        end
        
        mac = Expr(:macrocall,symbol("@gendict"),varname,:Uncertain,idxsets...)
        code = quote 
            $mac
            $code
            nothing
        end
        return code
    end
end

