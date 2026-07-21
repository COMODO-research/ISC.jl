module ISC

using Comodo
using Comodo.GeometryBasics
using Comodo.LinearAlgebra
using Printf
using DelimitedFiles

# Export imported modules for later possible use
export Comodo

# Export functions
include("ISC_misc.jl")
include("ISC_Abaqus.jl")
include("ISC_Checks.jl")
include("ISC_Geometry.jl")
include("ISC_Component.jl")
include("ISC_Optimization.jl")

# Misc functions
export iscdir

# export component functions
export comp_bolt, comp_WebStiffener, comp_SidePlate, comp_FL_v1, comp_IBeam_v1, comp_FL_Pantelis_v1

end #module ISC
