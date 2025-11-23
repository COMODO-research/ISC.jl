using ISC 
using Comodo
using Comodo.GLMakie

### INPUT
r           =       1.0                     # Corner Radii
g           =       1.0                     # Tolerance for Teeth
Gap         =       1.0                     # Gap for FL - FL
TH          =       10.0                    # Tooth Height      # at max    18      min ---
TD          =       11.5                    # Tooth Depth       # at max    20      min ---
TG          =       11.0                    # Tooth Gap         # at max    unl     min 11  
bFL         =       50.0                    # Width of FL       # at max    unl     min 34 
tFL         =       6.0                     # FL Thickness
nt_FL       =       3                       # no. of Teeth each side
L           =       98.0                    # Length of the Part
searchTol   =       1e-6  

Eh,Vh       =       mesh_flange(; r =1.0, g =1.0 ,Gap = 1.0  ,TH = 10.0   ,TD =  11.5,TG =  11.0,bFL = 50.0  ,tFL = 6.0, nt_FL =  3 ,L =  98.0 ,searchTol=  1e-6)    
#########################################    ## Visualization
Fh          =       element2faces(Eh)
markersize  = 25
linewidth   = 1
fig = Figure(size = (1200, 800))
# ax1 = Axis(fig[1, 1], xlabel = "X", ylabel = "Y")
ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z")
# scatter!(ax1, V1     , markersize = markersize *1   , color = :Green)
poly!(ax1, GeometryBasics.Mesh(Vh, Fh), strokewidth = linewidth, color = :white, strokecolor = :black, shading = FastShading, transparency = false)
# normalplot(ax1,Fh,Vh)
# wireframe!(ax1, GeometryBasics.Mesh(Vh, Eh), linewidth = 3, color = :red)
fig