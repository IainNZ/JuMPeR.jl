#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2016: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# src/solve_msgs.jl
# Warning messages for solving RobustModels.
# Included by src/solve.jl
#-----------------------------------------------------------------------

warn_ip_error() = Base.warn("""
Problem while solving a JuMPeR RobustModel:
The mixed-integer problem couldn't be solved. This may be due to the
solver not supporting lazy constraint callbacks, which are added by
default. You can disable all cutting planes by passing the option
disable_cuts=true to solve.
Error message:""")

warn_ip_infeasible() = Base.warn("""
Problem while solving a JuMPeR RobustModel:
The mixed-integer problem is infeasible.""")

warn_ip_unbounded() = Base.warn("""
Problem while solving a JuMPeR RobustModel:
The mixed-integer problem is unbounded. This may be due to:
* The problem requiring additional cutting planes for it to become
  bounded. Consider adding bounds to variables.
* The uncertainty set is unexpectedly empty. If using the default
  BasicUncertaintySet, you can enable cutting planes with the
  prefer_cuts=true option, which will usually detect empty sets.
* The problem is actually unbounded. You can suppress this warning
  message with the option suppress_warnings=true.""")

warn_lp_infeasible() = Base.warn("""
Problem while solving a JuMPeR RobustModel:
The continuous problem is infeasible.""")

warn_lp_unbounded_noray() = Base.warn("""
Problem while solving a JuMPeR RobustModel:
The continuous problem is unbounded and no ray is available. Lack
of a ray is normally due to the particular solver selected, but
the unboundedness itself may be due to:
* The problem requiring additional cutting planes for it to become
  bounded. Consider adding bounds to variables.
* The uncertainty set is unexpectedly empty. If using the default
  BasicUncertaintySet, you can enable cutting planes with the
  prefer_cuts=true option, which will usually detect empty sets.
* The problem is actually unbounded. You can suppress this warning
  message with the option suppress_warnings=true.""")

warn_lp_unbounded_sameray() = Base.warn("""
Problem while solving a JuMPeR RobustModel:
The continuous problem is unbounded and a ray is available. An
attempt was made to add cuts using this ray as the value of the
decision variables but the same ray was returned by solver again,
so the robust problem is assumed to be unbounded. This may be due to:
* The problem requiring additional cutting planes for it to become
  bounded. Consider adding bounds to variables.
* The uncertainty set is unexpectedly empty. If using the default
  BasicUncertaintySet, you can enable cutting planes with the
  prefer_cuts=true option, which will usually detect empty sets.
* The problem is actually unbounded. You can suppress this warning
  message with the option suppress_warnings=true.""")

warn_lp_unbounded() = Base.warn("""
Problem while solving a JuMPeR RobustModel:
The continuous problem is unbounded. This may be due to:
* The problem requiring additional cutting planes for it to become
  bounded. Consider adding bounds to variables.
* The uncertainty set is unexpectedly empty. If using the default
  BasicUncertaintySet, you can enable cutting planes with the
  prefer_cuts=true option, which will usually detect empty sets.
* The problem is actually unbounded. You can suppress this warning
  message with the option suppress_warnings=true.""")
