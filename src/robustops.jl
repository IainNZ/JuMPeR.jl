#############################################################################
# JuMPeR
# Julia for Mathematical Programming - extension for Robust Optimization
# See http://github.com/IainNZ/JuMPeR.jl
#############################################################################
# robustops.jl
# All the overloads for the new robust types introduced by JuMPeR. We share
# the same loose ordering as JuMP, just extended with the three new types:
# 1. Number
# 2. Variable
# 3. AffExpr
# 4. QuadExpr <- DISREGARD
# 5. Uncertain
# 6. UAffExpr
# 7. FullAffExpr
#############################################################################

# Number
# Number--Uncertain
(+)(lhs::Number, rhs::Uncertain) = UAffExpr([rhs],[+1.],convert(Float64,lhs))
(-)(lhs::Number, rhs::Uncertain) = UAffExpr([rhs],[-1.],convert(Float64,lhs))
(*)(lhs::Number, rhs::Uncertain) = UAffExpr([rhs],[convert(Float64,lhs)], 0.)
(/)(lhs::Number, rhs::Uncertain) = error("Cannot divide by an uncertain")
# Number--UAffExpr        - handled by JuMP
# Number--FullAffExpr     - handled by JuMP

# Variable
# Variable--Uncertain
(+)(lhs::Variable, rhs::Uncertain) = FullAffExpr([lhs],[UAffExpr(1.)],UAffExpr(rhs,+1.))
(-)(lhs::Variable, rhs::Uncertain) = FullAffExpr([lhs],[UAffExpr(1.)],UAffExpr(rhs,-1.))
(*)(lhs::Variable, rhs::Uncertain) = FullAffExpr([lhs],[UAffExpr(rhs,1.0)],UAffExpr())
(/)(lhs::Variable, rhs::Uncertain) = error("Cannot divide a variable by an uncertain")
# Variable--UAffExpr
(+)(lhs::Variable, rhs::UAffExpr) = FullAffExpr([lhs],[UAffExpr(1.)],    rhs)
(-)(lhs::Variable, rhs::UAffExpr) = FullAffExpr([lhs],[UAffExpr(1.)],0.0-rhs)
(*)(lhs::Variable, rhs::UAffExpr) = FullAffExpr([lhs],[rhs],UAffExpr())
(/)(lhs::Variable, rhs::UAffExpr) = error("Cannot divide a variable by an expression")
# Variable--FullAffExpr
(+)(lhs::Variable, rhs::FullAffExpr) = FullAffExpr(vcat(rhs.vars,lhs),vcat(rhs.coeffs,UAffExpr(1.)), rhs.constant)
(-)(lhs::Variable, rhs::FullAffExpr) = FullAffExpr(vcat(rhs.vars,lhs),vcat([0.0-rhs.coeffs[i] for i=1:length(rhs.coeffs)],UAffExpr(1.)),0.0-rhs.constant)
(*)(lhs::Variable, rhs::FullAffExpr) = error("Cannot multiply variable and expression")
(/)(lhs::Variable, rhs::FullAffExpr) = error("Cannot divide variable by expression")

# AffExpr
# AffExpr--Uncertain
(+)(lhs::AffExpr, rhs::Uncertain) = FullAffExpr(lhs.vars,map(UAffExpr,lhs.coeffs), UAffExpr([rhs],[+1.],lhs.constant))
(-)(lhs::AffExpr, rhs::Uncertain) = FullAffExpr(lhs.vars,map(UAffExpr,lhs.coeffs), UAffExpr([rhs],[-1.],lhs.constant))
(*)(lhs::AffExpr, rhs::Uncertain) = FullAffExpr(lhs.vars,[UAffExpr(rhs,lhs.coeffs[i]) for i=1:length(lhs.vars)], UAffExpr(rhs,lhs.constant))
(/)(lhs::AffExpr, rhs::Uncertain) = error("Cannot divide affine expression by an uncertain")
# AffExpr-UAffExpr
(+)(lhs::AffExpr, rhs::UAffExpr) = FullAffExpr(copy(lhs.vars),map(UAffExpr,lhs.coeffs),lhs.constant+rhs)
(-)(lhs::AffExpr, rhs::UAffExpr) = FullAffExpr(copy(lhs.vars),map(UAffExpr,lhs.coeffs),lhs.constant-rhs)
(*)(lhs::AffExpr, rhs::UAffExpr) = FullAffExpr(copy(lhs.vars),[lhs.coeffs[i]*rhs for i=1:length(lhs.vars)],lhs.constant*rhs)
(/)(lhs::AffExpr, rhs::UAffExpr) = error("Cannot divide affine expression by an uncertain expression")
# AffExpr-FullAffExpr
(+)(lhs::AffExpr, rhs::FullAffExpr) = FullAffExpr(
  vcat(lhs.vars, rhs.vars),
  vcat(map(UAffExpr,lhs.coeffs), rhs.coeffs),
  lhs.constant + rhs.constant)
(-)(lhs::AffExpr, rhs::FullAffExpr) = FullAffExpr(
  vcat(lhs.vars, rhs.vars),
  vcat(map(UAffExpr,lhs.coeffs), [0.0-rhs.coeffs[i] for i=1:length(rhs.coeffs)]),
  lhs.constant - rhs.constant)
(*)(lhs::AffExpr, rhs::FullAffExpr) = error("Cannot multiply expressions")
(/)(lhs::AffExpr, rhs::FullAffExpr) = error("Cannot divide expressions")

# Uncertain
(-)(lhs::Uncertain) = UAffExpr([lhs],[-1.],0.)
# Uncertain--Number
(+)(lhs::Uncertain, rhs::Number) = (+)(   +rhs, lhs)
(-)(lhs::Uncertain, rhs::Number) = (+)(   -rhs, lhs)
(*)(lhs::Uncertain, rhs::Number) = (*)(    rhs, lhs)
(/)(lhs::Uncertain, rhs::Number) = (*)(1.0/rhs, lhs)
# Uncertain--Variable
(+)(lhs::Uncertain, rhs::Variable) = (+)(rhs, lhs)
(-)(lhs::Uncertain, rhs::Variable) = FullAffExpr([rhs],[UAffExpr(-1.)],UAffExpr(lhs,+1.))
(*)(lhs::Uncertain, rhs::Variable) = (*)(rhs, lhs)
(/)(lhs::Uncertain, rhs::Variable) = error("Cannot divide uncertain by variable")
# Uncertain--AffExpr
(+)(lhs::Uncertain, rhs::AffExpr) = (+)(rhs, lhs)
(-)(lhs::Uncertain, rhs::AffExpr) = FullAffExpr(rhs.vars, map(UAffExpr,-rhs.coeffs), UAffExpr([lhs],[1.],-rhs.constant))
(*)(lhs::Uncertain, rhs::AffExpr) = (*)(rhs, lhs)
(/)(lhs::Uncertain, rhs::AffExpr) = error("Cannot divide uncertain by expression")
# Uncertain--Uncertain
(+)(lhs::Uncertain, rhs::Uncertain) = UAffExpr(vcat(lhs,rhs),[1.,+1.],0.)
(-)(lhs::Uncertain, rhs::Uncertain) = UAffExpr(vcat(lhs,rhs),[1.,-1.],0.)
(*)(lhs::Uncertain, rhs::Uncertain) = error("Cannot multiply two uncertains")
(/)(lhs::Uncertain, rhs::Uncertain) = error("Cannot divide two uncertains")
# Uncertain--UAffExpr
(+)(lhs::Uncertain, rhs::UAffExpr) = UAffExpr(vcat(lhs,rhs.vars),vcat(1.0, rhs.coeffs), rhs.constant)
(-)(lhs::Uncertain, rhs::UAffExpr) = UAffExpr(vcat(lhs,rhs.vars),vcat(1.0,-rhs.coeffs),-rhs.constant)
(*)(lhs::Uncertain, rhs::UAffExpr) = error("Cannot multiply uncertain and expression")
(/)(lhs::Uncertain, rhs::UAffExpr) = error("Cannot divide uncertain by expression")
# Uncertain--FullAffExpr
(+)(lhs::Uncertain, rhs::FullAffExpr) = FullAffExpr(copy(rhs.vars),copy(rhs.coeffs),lhs+rhs.constant)
(-)(lhs::Uncertain, rhs::FullAffExpr) = FullAffExpr(copy(rhs.vars),[0.0-rhs.coeffs[i] for i=1:length(rhs.coeffs)],lhs-rhs.constant)
(*)(lhs::Uncertain, rhs::FullAffExpr) = error("Cannot multiply uncertainty by uncertain expression")
(/)(lhs::Uncertain, rhs::FullAffExpr) = error("Cannot divide uncertainty by uncertain expression")

# UAffExpr
# UAffExpr--Number        - handled by JuMP
# UAffExpr--Variable
(+)(lhs::UAffExpr, rhs::Variable) = (+)(rhs,lhs)
(-)(lhs::UAffExpr, rhs::Variable) = FullAffExpr([rhs],[UAffExpr(-1.)],lhs)
(*)(lhs::UAffExpr, rhs::Variable) = (*)(rhs,lhs)
(/)(lhs::UAffExpr, rhs::Variable) = error("Cannot divide by variable")
# UAffExpr--AffExpr
(+)(lhs::UAffExpr, rhs::AffExpr) = (+)(rhs,lhs)
(-)(lhs::UAffExpr, rhs::AffExpr) = (+)(0.0-rhs,lhs)
(*)(lhs::UAffExpr, rhs::AffExpr) = (*)(rhs,lhs)
(/)(lhs::UAffExpr, rhs::AffExpr) = error("Cannot divide by affine expression")
# UAffExpr--Uncertain
(+)(lhs::UAffExpr, rhs::Uncertain) = (+)(rhs,lhs)
(-)(lhs::UAffExpr, rhs::Uncertain) = UAffExpr(vcat(rhs,lhs.vars),vcat(-1.0,lhs.coeffs),lhs.constant)
(*)(lhs::UAffExpr, rhs::Uncertain) = (*)(rhs,lhs)
(/)(lhs::UAffExpr, rhs::Uncertain) = error("Cannot divide by uncertain")
# UAffExpr--UAffExpr
(+)(lhs::UAffExpr, rhs::UAffExpr) = UAffExpr(vcat(lhs.vars,rhs.vars),vcat(lhs.coeffs, rhs.coeffs),lhs.constant+rhs.constant)
(-)(lhs::UAffExpr, rhs::UAffExpr) = UAffExpr(vcat(lhs.vars,rhs.vars),vcat(lhs.coeffs,-rhs.coeffs),lhs.constant-rhs.constant)
(*)(lhs::UAffExpr, rhs::UAffExpr) = error("Cannot multiply two expressions")
(/)(lhs::UAffExpr, rhs::UAffExpr) = error("Cannot divide two expressions")
# UAffExpr--FullAffExpr
(+)(lhs::UAffExpr, rhs::FullAffExpr) = FullAffExpr(rhs.vars,rhs.coeffs,lhs+rhs.constant)
(-)(lhs::UAffExpr, rhs::FullAffExpr) = FullAffExpr(rhs.vars,[0.0-c for c in rhs.coeffs],lhs-rhs.constant)
(*)(lhs::UAffExpr, rhs::FullAffExpr) = error("Cannot multiply two expressions")
(/)(lhs::UAffExpr, rhs::FullAffExpr) = error("Cannot divide two expressions")

# FullAffExpr
# FullAffExpr--Number     - handled by JuMP
# FullAffExpr--Variable
(+)(lhs::FullAffExpr, rhs::Variable) = (+)(rhs,lhs)
(-)(lhs::FullAffExpr, rhs::Variable) = FullAffExpr(vcat(lhs.vars,rhs),vcat(lhs.coeffs,UAffExpr(-1.)), lhs.constant)
(*)(lhs::FullAffExpr, rhs::Variable) = error("Cannot")
(/)(lhs::FullAffExpr, rhs::Variable) = error("Cannot")
# FullAffExpr--AffExpr
(+)(lhs::FullAffExpr, rhs::AffExpr) = (+)(rhs,lhs)
(-)(lhs::FullAffExpr, rhs::AffExpr) = FullAffExpr(
  vcat(lhs.vars, rhs.vars),
  vcat(lhs.coeffs, -rhs.coeffs),
  lhs.constant - rhs.constant)
(*)(lhs::FullAffExpr, rhs::AffExpr) = error("Cannot")
(/)(lhs::FullAffExpr, rhs::AffExpr) = error("Cannot")
# FullAffExpr--Uncertain
(+)(lhs::FullAffExpr, rhs::Uncertain) = (+)(rhs,lhs)
(-)(lhs::FullAffExpr, rhs::Uncertain) = FullAffExpr(copy(lhs.vars),copy(lhs.coeffs),lhs.constant-rhs)
(*)(lhs::FullAffExpr, rhs::Uncertain) = error("Cannot")
(/)(lhs::FullAffExpr, rhs::Uncertain) = error("Cannot")
# FullAffExpr--UAffExpr
(+)(lhs::FullAffExpr, rhs::UAffExpr) = (+)(rhs,lhs)
(-)(lhs::FullAffExpr, rhs::UAffExpr) = FullAffExpr(lhs.vars,lhs.coeffs,lhs.constant-rhs)
(*)(lhs::FullAffExpr, rhs::UAffExpr) = error("Cannot")
(/)(lhs::FullAffExpr, rhs::UAffExpr) = error("Cannot")
# FullAffExpr--FullAffExpr
(+)(lhs::FullAffExpr, rhs::FullAffExpr) = FullAffExpr(vcat(lhs.vars,rhs.vars),vcat(lhs.coeffs,rhs.coeffs),lhs.constant+rhs.constant)
(-)(lhs::FullAffExpr, rhs::FullAffExpr) = FullAffExpr(vcat(lhs.vars,rhs.vars),vcat(lhs.coeffs,[0.0-c for c in rhs.coeffs]),lhs.constant-rhs.constant)
(*)(lhs::FullAffExpr, rhs::FullAffExpr) = error("Cannot")
(/)(lhs::FullAffExpr, rhs::FullAffExpr) = error("Cannot")

# Constraints
# Number, Variable, AffExpr
for (sgn, osgn) in ( (:<=,:>=), (:(==),:(==)), (:>=,:<=) )
    for old_typ in (:Number, :Variable, :AffExpr)
        for new_typ in (:Uncertain, :UAffExpr, :FullAffExpr)
            @eval $(sgn)(lhs::$(old_typ), rhs::$(new_typ)) = $(osgn)(rhs, lhs)
        end
    end
end
# Uncertain, UAffExpr, FullAffExpr vs Number
(<=)(lhs::Uncertain, rhs::Number)   = UncSetConstraint(UAffExpr(lhs,1.0), -Inf,  rhs)
(==)(lhs::Uncertain, rhs::Number)   = UncSetConstraint(UAffExpr(lhs,1.0),  rhs,  rhs)
(>=)(lhs::Uncertain, rhs::Number)   = UncSetConstraint(UAffExpr(lhs,1.0),  rhs, +Inf)
(<=)(lhs::UAffExpr, rhs::Number)    = UncSetConstraint(lhs, -Inf, rhs - lhs.constant)
(==)(lhs::UAffExpr, rhs::Number)    = UncSetConstraint(lhs,  rhs - lhs.constant, rhs - lhs.constant)
(>=)(lhs::UAffExpr, rhs::Number)    = UncSetConstraint(lhs, rhs - lhs.constant, +Inf)
(<=)(lhs::FullAffExpr, rhs::Number) = UncConstraint(lhs, -Inf,  rhs - lhs.constant.constant)
(==)(lhs::FullAffExpr, rhs::Number) = UncConstraint(lhs,  rhs - lhs.constant.constant,  rhs - lhs.constant.constant)
(>=)(lhs::FullAffExpr, rhs::Number) = UncConstraint(lhs,  rhs - lhs.constant.constant, +Inf)
# Uncertain, UAffExpr, FullAffExpr vs Other
for sgn in ( :<=, :(==), :>= )
    for new_typ in (:Uncertain, :UAffExpr, :FullAffExpr)
        for other_typ in (:Variable, :AffExpr, :Uncertain, :UAffExpr, :FullAffExpr)
            @eval $(sgn)(lhs::$(new_typ), rhs::$(other_typ)) = $(sgn)(lhs - rhs, 0.0)
        end
    end
end

#############################################################################
# High-level operators
# Currently supported
#  - sum
#  - dot
#############################################################################

# SUM
# Neither is strictly needed, but more efficient than fallbacks in JuMP
Base.sum(j::JuMPArray{Uncertain}) = UAffExpr(vec(j.innerArray), ones(length(j.innerArray)), 0.0)
Base.sum(j::JuMPDict{Uncertain})  = UAffExpr(collect(values(j.tupledict)), ones(length(j.tupledict)), 0.0)

# DOT
if VERSION < v"0.4-"
    import JuMP.vecdot
    export vecdot
end