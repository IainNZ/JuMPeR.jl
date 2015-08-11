# This file attempts to use as much of the machinery provided by JuMP
# to provide the @defUnc functionality for Uncertains, just like @defVar
# does for variables. @defUnc is essentially a direct copy of @defVar,
# with the differences clearly marked. If something is not marked,
# assume it is a direct copy.
using Base.Meta

macro defUnc(args...)
    length(args) <= 1 &&
        error("in @defUnc: expected model as first argument, then uncertain parameter information.")
    m = args[1]
    x = args[2]
    extra = vcat(args[3:end]...)
    m = esc(m)

    t = :Cont
    gottype = 0
    # Identify the variable bounds. Five (legal) possibilities are "x >= lb",
    # "x <= ub", "lb <= x <= ub", "x == val", or just plain "x"
    if isexpr(x,:comparison)
        # We have some bounds
        if x.args[2] == :>= || x.args[2] == :≥
            if length(x.args) == 5
                # ub >= x >= lb
                x.args[4] == :>= || x.args[4] == :≥ || error("Invalid uncertain parameter bounds")
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
                    error("in @defUnc ($var): expected <= operator after uncertain parameter name.")
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
        elseif x.args[2] == :(==)
            # fixed variable
            var = x.args[1]
            @assert length(x.args) == 3
            lb = esc(x.args[3])
            ub = esc(x.args[3])
            gottype = 1
            t = :Fixed
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

    # separate out keyword arguments
    kwargs = filter(ex->isexpr(ex,:kw), extra)
    extra = filter(ex->!isexpr(ex,:kw), extra)

    # process keyword arguments
    #------------------------------------------------------------------
    # MODIFIED TO DISABLE COLUMN GENERATION
    for ex in kwargs
        error("in @defUnc ($unc): Unrecognized keyword argument $(ex.args[1])")
    end
    #------------------------------------------------------------------

    sdp = any(t -> (t == :SDP), extra)
    symmetric = sdp | any(t -> (t == :Symmetric), extra)
    extra = filter(x -> (x != :SDP && x != :Symmetric), extra) # filter out SDP and sym tag
    #------------------------------------------------------------------
    # DON'T ALLOW SDP/SYMMETRIC
    if sdp || symmetric
        error("in @defUnc ($var): SDP and Symmetric not supported.")
    end
    #------------------------------------------------------------------

    # Determine variable type (if present).
    # Types: default is continuous (reals)
    if length(extra) > 0
        if t == :Fixed
            error("in @defUnc ($var): unexpected extra arguments when declaring a fixed uncertain parameter.")
        end
        if extra[1] in [:Bin, :Int, :SemiCont, :SemiInt]
            gottype = 1
            t = extra[1]
        end

        if t == :Bin
            if (lb != -Inf || ub != Inf) && !(lb == 0.0 && ub == 1.0)
            error("in @defUnc ($var): bounds other than [0, 1] may not be specified for binary uncertain parameters.\nThese are always taken to have a lower bound of 0 and upper bound of 1.")
            else
                lb = 0.0
                ub = 1.0
            end
        end

        gottype == 0 &&
            error("in @defUnc ($var): syntax error")
    end

    # Handle the column generation functionality
    #------------------------------------------------------------------
    # SKIP
    #------------------------------------------------------------------

    if isa(var,Symbol)
        # Easy case - a single variable
        #------------------------------------------------------------------
        # CANNOT OCCUR
        # sdp && error("Cannot add a semidefinite scalar variable")
        #------------------------------------------------------------------
        return assert_validmodel(m, quote
            # CHANGED TO UNCERTAIN
            $(esc(var)) = Uncertain($m,$lb,$ub,$(quot(t)),$(string(var)))
            # TODO
            # registervar($m, $(quot(var)), $(esc(var)))
        end)
    end
    isa(var,Expr) || error("in @defUnc: expected $var to be an uncertain parameter name")


    # We now build the code to generate the variables (and possibly the JuMPDict
    # to contain them)
    refcall, idxvars, idxsets, idxpairs, condition = buildrefsets(var)
    #------------------------------------------------------------------
    # IGNORE SYMMETRIC CODE, VARIABLE -> UNCERTAIN
    code = :( $(refcall) = Uncertain($m, $lb, $ub, $(quot(t)), "") )
    looped = getloopedcode(var, code, condition, idxvars, idxsets, idxpairs, :Uncertain)
    varname = esc(getname(var))
    return assert_validmodel(m, quote
        $looped
        push!($(m).ext[:Robust].dictList, $varname)
        #registervar($m, $(quot(getname(var))), $varname)
        $varname
    end)
    #------------------------------------------------------------------
end


function JuMP.constructconstraint!(faff::FullAffExpr, sense::Symbol)
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