module OneSiteZ2AD
using LinearAlgebra, Random, Zygote, JSON, SHA
module Legacy
include("../split/network.jl")
include("../split/objective.jl")
include("../split/optimize.jl")
end
include("parameters.jl")
include("state_io.jl")
include("objective.jl")
include("optimize.jl")
include("runner.jl")
end
