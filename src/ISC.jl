module ISC

using Comodo
using Comodo.GeometryBasics

# Export imported modules for later possible use
# export Comodo

# Export functions
include("functions.jl")
export iscdir, mesh_sideplate

end # module ISC
