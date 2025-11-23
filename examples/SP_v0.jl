using ISC 
using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Statistics
using Comodo.Rotations
using Comodo.LinearAlgebra
using Printf

## INPUTS
B_SP    = 20.0
TG      = 11.0        # to be linked
TH      = 15.0        # to be linked
TThc     = 6.0 
g       = 0.5       # to be linked #gap on each
tSP     = 6.0    
mesh_spacing   = 4.5
nSections      = 2                    # no. of teeth in each flange

searchTol = 1e-6

EN6, VN6 = mesh_sideplate(; B_SP=B_SP, TG=TG, TH=TH, g=g, TThc=TThc, tSP=tSP, mesh_spacing=mesh_spacing, nSections=nSections, searchTol = searchTol)    

FN6       = element2faces(EN6)
FN6_Bound = boundaryfaces(FN6)

nSectionsHoles = 2*nSections
Ht_SP = nSectionsHoles*(TG + TH) + TG -2*g

################################################  # make Set-3 for the BC SP Corners VN6 
V_Set_3 = filter(p ->   (isapprox(p[1], 0 + (0.5 * B_SP), atol=searchTol) || 
                        isapprox(p[1], 0 - (0.5 * B_SP), atol=searchTol)) &&
                        (isapprox(p[2], 0 + (0.5 * TH + 0.5 * TG + 0.5 * TG - g), atol=searchTol) || 
                        isapprox(p[2], 0 - (Ht_SP - (0.5 * TH + 0.5 * TG + 0.5 * TG - g)), atol=searchTol)) &&
                        (isapprox(p[3], 0, atol=searchTol) || 
                        isapprox(p[3], tSP, atol=searchTol)), VN6)

node_num_V_Set_3 = findall(p -> p in V_Set_3, VN6)
##################################################    ## Creating Input File


file2       =       joinpath(iscdir(), "assets", "temp", "temp_SP.inp")
open(file2,"w") do io
    println(io,"*Part, name=SP")
    println(io,"* Node")
    for (i,v) in enumerate(VN6)
        println(io,@sprintf("%i", i) * ", " * join([@sprintf("%.16e",x) for x in v],", "))
    end

    println(io,"*Element, type=C3D8R")
    for (i,e) in enumerate(EN6)
        println(io,@sprintf("%i", i) * ", " * join([@sprintf("%i",x) for x in e],", "))
    end
end
##   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   
# Visualization

FN6p, VN6p = separate_vertices(FN6, VN6)

markersize = 25
linewidth  = 1
fig = Figure(size = (800, 1200))
ax1 = AxisGeom(fig[1, 1], title="Side plate design and mesh")

# Fp = F1
# Vp = V1
# Ep = boundaryedges(Fp)
# poly!(ax1, GeometryBasics.Mesh(Vp, Fp), strokewidth = linewidth, color = :white, strokecolor = :black, shading = FastShading, transparency = false)
# normalplot(ax1,Fp,Vp)
# wireframe!(ax1, GeometryBasics.Mesh(Vp, Ep), linewidth = 10, color = :red)

meshplot!(ax1, FN6p, VN6p, strokewidth = linewidth)

# scatter!(ax1, V1     , markersize = markersize *1   , color = :Green)
# poly!(ax1, GeometryBasics.Mesh(VN6, FN6), strokewidth = linewidth, color = :white, strokecolor = :black, shading = FastShading, transparency = false)
# normalplot(ax1,FN6,VN6)
# wireframe!(ax1,GeometryBasics.Mesh(VN5,boundaryedges(FN5) ), linewidth=5,color=:red)

fig





