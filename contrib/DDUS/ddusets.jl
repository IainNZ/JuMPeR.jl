######
# Data-Driven Sets for JuMPeR
######

module DDUSets

using JuMPeR  #VG Talk to Iain about removing this

include("helpers.jl")
include("FBOracle.jl")
include("UIOracle.jl")
include("UCSOracle.jl")
include("LCXOracle.jl")
include("UMOracle.jl")
#include("UDYOracle.jl")  VG Needs SDP support...  add back later
end