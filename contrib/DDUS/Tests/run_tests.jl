###
# Test harness for the DDUSets tests
###

include("../ddusets.jl")  #To add the module to the path.  Talk to Iain to remove

using DDUSets
using FactCheck
using JuMPeR

include("test_helpers.jl")  #loads up generic functionality

include("StatHelpers_tests.jl")
include("FBOracle_tests.jl")
include("UIOracle_tests.jl")
include("UCSOracle_tests.jl")  #retains a dependence on Unsetsolver
include("LCXOracle_tests.jl")
include("UMOracle_tests.jl")
#include("UDYOracle_tests.jl")  #needs to be made compatible toFact Check