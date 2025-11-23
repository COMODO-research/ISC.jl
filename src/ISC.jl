module ISC

using Comodo
using Comodo.GeometryBasics

# Export imported modules for later possible use
# export Comodo

# Export functions
include("functions.jl")
export iscdir, mesh_sideplate, mesh_flange, f1, f2

end #module ISC
