#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# Wranglers
# Wranglers are "robustifying operators". They take a list of constraints,
# some of which maybe null or partially defined, and are reponsible for
# "robustifying" them by either reformulating them, providing a cutting plane
# algorithm, and/or a sampling procedure. They have the full statement of
# the uncertain problem available to them, in particular the uncertainty set
# and uncertain bounds.
#############################################################################

#############################################################################
# AbstractWrangler
# All wranglers implement the interface defined by AbstractWrangler
abstract AbstractWrangler

# registerConstraint
# Notifies the wrangler that it is responsible for this constraint, and 
# passes any preferences provided via the solve command. Returns a dictionary
# where the keys are the symbols :Cut, :Reform, and :Sample and the values
# are true or false, where true indicates the wrangler has selected these
# operations for this constraint. (not currently used)
registerConstraint(w::AbstractWrangler, con, ind::Int, prefs) = error("Not implemented!")

# setup
# Gives wrangler time to do any setup it needs to do. Called after all
# constraints have been registered. Examples of work that could be done here
# include transforming the uncertainty set and generating a cutting plane
# model. May be called multiple times - this should be handled by the wrangler
setup(w::AbstractWrangler, rm::Model) = error("Not implemented")

# generateCut
# Called in the main loop every iteration. m is the actual current model, aka
# the master model, that will have the current solution and to which 
# constraints should be added. Returns number of cuts added.
generateCut(w::AbstractWrangler, rm::Model, ind::Int, m::Model) = error("Not implemented")

# generateReform
# Called before the main loop, adds anything it wants to the model
generateReform(w::AbstractWrangler, rm::Model, ind::Int, m::Model) = error("Not implemented")

#############################################################################
# Default included wranglers
include("wrangler_poly.jl")         # PolyhedralWrangler
include("wrangler_bertsim.jl")      # BertSimWrangler