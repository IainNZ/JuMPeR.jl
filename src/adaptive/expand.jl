#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2015: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# src/adaptive/expand.jl
# Adaptive robust optimization support - pre-solve expansion
#-----------------------------------------------------------------------

any_adaptive(u::UncAffExpr) = any(v->isa(v,Adaptive), u.vars)

function expand_adaptive(rm::Model)
    rd = getRobust(rm)::RobustData

    length(rd.adpPolicy) == 0 && return

    new_cons = UncConstraint[]

    subs = Any[]
    for i in 1:rd.numAdapt
        pol = rd.adpPolicy[i]
        if pol == :Static
            # Static policy = no dependence
            x_static = Variable(rm,rd.adpLower[i],rd.adpUpper[i],rd.adpCat[i],rd.adpNames[i])
            push!(subs, x_static)
        elseif pol == :Affine
            deps = rd.adpDependsOn[i]
            if rd.adpStage[i] != 0                
                # Stage-wise dependence being used
                error("Not implemetned")
            else
                # Explicit dependency
                # Create auxiliary variables
                aff_pol = Dict()
                for j in eachindex(deps)
                    vname = utf8(string(adp_str(rm,i), "{", deps[j], "}"))
                    aff_pol[j] = Variable(rm,-Inf,+Inf,:Cont,vname)
                end
                aff_con = Variable(rm,-Inf,+Inf,:Cont,string(adp_str(rm,i), utf8("{}")))
                # Build the policy
                x_aff = UncAffExpr()
                for j in eachindex(deps)
                    push!(x_aff, UncExpr(deps[j]), aff_pol[j])
                end
                push!(x_aff, UncExpr(1), aff_con)
                push!(subs, x_aff)
                # Add bound constraints on the policy
                if rd.adpLower[i] != -Inf && rd.adpUpper[i] != +Inf
                    push!(new_cons, UncConstraint(x_aff, rd.adpLower[i], rd.adpUpper[i]))
                end
            end
        else
            error("errr")
        end
    end


    # Replace the other constraints
    init_n = length(rd.uncertainconstr)
    for i in 1:init_n
        uncaffcon = rd.uncertainconstr[i]
        lhs = uncaffcon.terms
        !any_adaptive(lhs) && continue
        new_lhs = UncAffExpr()
        for (coeff, var) in lhs
            if isa(var, Adaptive)
                new_var = subs[var.id]
                new_lhs += coeff * new_var
            else
                push!(new_lhs, coeff, var)
            end
        end
        new_lhs.constant = lhs.constant
        push!(new_cons, UncConstraint(new_lhs, uncaffcon.lb, uncaffcon.ub))
        lhs.vars = JuMPeRVar[]
        lhs.coeffs = UncExpr[]
        lhs.constant = UncExpr()
    end
    

    # Replace the number-and-adaptive constraints
    for varaffcon in rd.varaffcons
        lhs = varaffcon.terms
        new_lhs = UncAffExpr()
        for (coeff, var) in lhs
            if isa(var, Adaptive)
                new_var = subs[var.id]
                new_lhs += coeff * new_var
            else
                push!(new_lhs, coeff, var)
            end
        end
        push!(new_cons, UncConstraint(new_lhs, varaffcon.lb, varaffcon.ub))
    end

    map(c->addConstraint(rm, c), new_cons)

    #println(rm)
end