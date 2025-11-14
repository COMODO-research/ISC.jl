using ISC 
using Comodo
using Comodo.GLMakie

## INPUTS
B_SP    = 20.0
TG      = 11.0        # to be linked
TH      = 15.0        # to be linked
TThc     = 6.0 
g       = 0.5       # to be linked #gap on each
tSP     = 6.0    
mesh_spacing   = 4.5
nSections      = 2                    # no. of teeth in each flange

EN6, VN6 = mesh_sideplate(; B_SP=B_SP, TG=TG, TH=TH, g=g, TThc=TThc, tSP=tSP, mesh_spacing=mesh_spacing, nSections=nSections)    

# Visualization
FN6       = element2faces(EN6)
FN6_Bound = boundaryfaces(FN6)
FN6_Boundp, VN6p = separate_vertices(FN6_Bound, VN6)

markersize = 25
linewidth  = 1

n = length(EN6)

fig = Figure(size = (1200, 800))
ax1 = AxisGeom(fig[1, 1], title="Side plate design and mesh, featuring $n elements")
hp = meshplot!(ax1, FN6_Boundp, VN6p, strokewidth = linewidth)
fig