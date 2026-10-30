module ToricTransferSpectrum
include("../simple-update/one-site-c4-z2/core.jl")
using .OneSiteZ2SU
include("numerics.jl")
export build_transfer
include("input.jl")
include("../ad/one-site-c4-z2/core.jl")
include("output.jl")
end
