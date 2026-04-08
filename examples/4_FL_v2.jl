# using ISC
using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Statistics
using Comodo.Rotations
using Comodo.LinearAlgebra
using Printf
### Going to create the I Beam.
### INPUT
r           =       1.0                     # Corner Radii
g           =       1.0                     # Tolerance for Teeth
Gap         =       1.0                     # Gap for FL - FL
TH          =       10.0                    # Tooth Height      # at max    18      min ---
TD          =       11.5                    # Tooth Depth       # at max    20      min ---
TG          =       11.0                    # Tooth Gap         # at max    unl     min 11  
bFL         =       50.0                    # Width of FL       # at max    unl     min 34 
tFL         =       6.0                     # FL Thickness
depthWeb    =       50.0
thcWeb      =       6.0
nt_FL       =       3                       # no. of Teeth each side
L           =       98.0                    # Length of the Part
searchTol   =       1e-6 
### DERIVED INPUTS
shiftDist   =       (TH + TG)
FL_B_L      =       L-(nt_FL *(shiftDist) + 0.5*TG - 0.5*Gap)
beFL        =       bFL - 2 * TD            # Effective W.of FL
nSections   =       nt_FL
n           =       5                                       # n=ceil(Int,((0.5*TH+r)/1.5))
### ### ### FL_TOP STARTS
### 0th LOFTING - EL 0 - HI to ED           
e0_h        =       0.5*TH+r
e0_w        =       0.5*thcWeb
V0_S        =       Vector{Point{3, Float64}}(undef, n)
v_1         =       Point{3, Float64}(0.0, 0.0 , 0.0)
v_2         =       Point{3, Float64}(0.0, e0_h , 0.0)
for (i, w) in enumerate(range(0.0, 1.0, n))             
    V0_S[i] =       (1.0 - w) * v_1 + w * v_2          
end
V0_E        =       [Point{3, Float64}(v[1]-e0_w, v[2], v[3]) for v in V0_S]
F0, V0      =       loftlinear(V0_S, V0_E; num_steps = 10, close_loop = false, face_type = :quad)
### 1st LOFTING - EL 1 - HI to AB
e1_h        =       0.5*TH+r
e1_w        =       0.5*beFL-0.5*TH-e0_w  
V1_S        =       [Point{3, Float64}(v[1], v[2], v[3]) for v in V0_E]                               
V1_E        =       [Point{3, Float64}(v[1]-e1_w, v[2], v[3]) for v in V1_S]
F1, V1      =       loftlinear(V1_S, V1_E; num_steps = 10, close_loop = false, face_type = :quad)
### 2nd LOFTING - EL 2 - AB BC to XY       
e2_w        =       0.5*TH 
Vr          =       [Point{3, Float64}(r * cosd(t), r * sind(t), 0.0) for t in range(0.0, 90, n+(n-1))]
Vr          =       [Point{3, Float64}(v[1]-e0_w-e1_w-e1_h, v[2], v[3]) for v in Vr]
Vb          =       Vector{Point{3, Float64}}(undef, n)  
Vb1_E       =       [Point{3, Float64}(v[1], v[2], v[3]) for v in V1_E]
Vb2_E       =       Vector{Point{3, Float64}}(undef, n)          
v_1         =       Point{3, Float64}(-e0_w-e1_w        , e1_h , 0.0)
v_2         =       Point{3, Float64}(-e0_w-e1_w-e1_h   , e1_h , 0.0)
for (i, w) in enumerate(range(0.0, 1.0, n))              
    Vb2_E[i] = (1.0 - w) * v_1 + w * v_2                
end
Vb          =       [Vb1_E;Vb2_E]
Vb          =       unique(Vb, dims=1)
F2, V2      =       loftlinear(Vb, Vr; num_steps = 10, close_loop = false, face_type = :quad)
### 3rd LOFTING - EL 3 - HX to GF
e3_h        =       0.5*TG-r   
V3_S0       =       filter(p -> isapprox(p[2], 0.0 ,atol=searchTol), V0)
V3_S1       =       filter(p -> isapprox(p[2], 0.0 ,atol=searchTol), V1)
V3_S2       =       filter(p -> isapprox(p[2], 0.0 ,atol=searchTol), V2)
V3_S        =       [V3_S0;V3_S1;V3_S2]
V3_S        =       unique(V3_S, dims=1)
V3_E        =       [Point{3, Float64}(v[1], v[2]-e3_h, v[3]) for v in V3_S]
F3, V3      =       loftlinear(V3_S, V3_E; num_steps = 5, close_loop = false, face_type = :quad)
### 4th LOFTING - EL 4 - YC to DE
e4_w        =       TD-r
V4_S        =       filter(p -> isapprox(p[1], -e0_w-e1_w-e1_h ,atol=searchTol), V2)
V4_E        =       [Point{3, Float64}(v[1]-e4_w, v[2], v[3]) for v in V4_S]
F4, V4      =       loftlinear(V4_E, V4_S; num_steps = n, close_loop = false, face_type = :quad)
### 5th Merging 1st 2nd 3rd 4th
V5          =       [V0; V1; V2; V3; V4]
F5          =       [F0;
                    [f .+ length(V0) for f in F1];
                    [f .+ length(V0) .+ length(V1) for f in F2];
                    [f .+ length(V0) .+ length(V1) .+ length(V2) for f in F3]
                    [f .+ length(V0) .+ length(V1) .+ length(V2).+ length(V3) for f in F4]]
F5, V5      =       mergevertices(F5, V5)
Eb5         =       boundaryedges(F5)                  
V5        .+=       Point{3, Float64}(0.0, +e3_h, 0.0)
### 6th Copy of basic element i.e. V5
V6          =       [V5;
                    [Point{3, Float64}(-v[1],  v[2], v[3]) for v in V5];
                    [Point{3, Float64}( v[1], -v[2], v[3]) for v in V5];
                    [Point{3, Float64}(-v[1], -v[2], v[3]) for v in V5]]    
F6          =       [F5;
                    [reverse(f)  .+    length(V5)      for f in F5];
                    [reverse(f)  .+    length(V5)*2    for f in F5];
                    [f           .+    length(V5)*3    for f in F5]]
### 7th Copy of nSections
nf          =       length(F6)
nv          =       length(V6)
F7          =       Vector{QuadFace{Int64}}(undef,nf*nSections)
V7          =       Vector{Point{3,Float64}}(undef,nv*nSections)
i_f         =       1
i_v         =       1
s           =       0
for q in 1:nSections
    F7[i_f:i_f+nf-1]  = [f.+s for f in F6]
    V7[i_v: i_v+nv-1] = [Point{3, Float64}(v[1], v[2] - (q-1)*shiftDist, v[3]) for v in V6]
    global i_f += nf
    global i_v += nv
    global s   += nv
end
F7, V7      =      mergevertices(F7, V7)
### 8th Loft Top-Part_1 (Left)
V8          =       [[Point{3, Float64}(v[1],  -v[2]  .+ shiftDist-0.5*TG+r, v[3]) for v in V1];
                    [Point{3, Float64}(v[1],   -v[2]  .+ shiftDist-0.5*TG+r, v[3]) for v in V2];
                    [Point{3, Float64}(v[1],   -v[2]  .+ shiftDist-0.5*TG+r, v[3]) for v in V4]]

F8          =       [[reverse(f) for f in F1];
                    [reverse(f) .+ length(V1) for f in F2]; 
                    [reverse(f) .+ length(V1) .+ length(V2) for f in F4]]
F8, V8      =       mergevertices(F8, V8)
### 9th Loft Top-Part_2 (Left)
y           =       0.5*TG - 0.5*g - r
V9_S        =       filter(p -> isapprox(p[2], 0.5*TG + TH + r ,atol=searchTol), V8)
indSort     =       reverse(sortperm([v[1] for v in V9_S]))
V9_S        =       V9_S[indSort]
V9_E        =       [Point{3, Float64}(v[1],  v[2]  .+ y, v[3]) for v in V9_S]
indSort     =       reverse(sortperm([v[1] for v in V9_S]))
V9_E        =       V9_E[indSort]
F9,V9       =       loftlinear(V9_E, V9_S ; num_steps = 10, close_loop = false, face_type = :quad)
### 10th Loft Top-Part_3 (Center Left)
V10_S1      =       filter(p -> isapprox(p[1], -e0_w ,atol=searchTol), V8)
V10_S2      =       filter(p -> isapprox(p[1], -e0_w ,atol=searchTol), V9)
V10_S       =       [V10_S1;V10_S2]
indSort     =       reverse(sortperm([v[2] for v in V10_S]))
V10_S       =       V10_S[indSort]
V10_S       =       unique(V10_S, dims=1)
V10_E       =       [Point{3, Float64}(v[1]+e0_w, v[2], v[3]) for v in V10_S]
indSort     =       reverse(sortperm([v[2] for v in V10_E]))
V10_E       =       V10_E[indSort]
F10,V10     =       loftlinear(V10_S, V10_E; num_steps = 10, close_loop = false, face_type = :quad)
### 11th Top-Part Merge 8th, 9th and 10th
V11         =       [V8; V9; V10]                    
F11         =       [F8;
                    [f .+ length(V8) for f in F9]
                    [f .+ length(V8) .+ length(V9) for f in F10]]
F11, V11    =  mergevertices(F11, V11)
Eb11        =       boundaryedges(F11)        #??
### Loft Top-Part Merge & Mirror 
V_N2        =       [V11;
                    [Point{3, Float64}(-v[1],  v[2], v[3]) for v in V11]]     # top right    # bottom left
F_N2        =       [F11;
                    [reverse(f)  .+    length(V11)     for f in F11]]
F_N2, V_N2  =       mergevertices(F_N2, V_N2)
Eb_N2       =       boundaryedges(F_N2)
### Loft Bottom-Part
V_N3_S      =       filter(p -> isapprox(p[2],-(nSections-0.5)*(TG+TH),atol=searchTol), V7)
indSort     =       reverse(sortperm([v[1] for v in V_N3_S]))
V_N3_S      =       V_N3_S[indSort]
V_N3_E      =       [Point{3, Float64}(v[1], -(FL_B_L-0.5*TH)+v[2] , v[3]) for v in V_N3_S]
F_N3, V_N3  =       loftlinear(V_N3_S, V_N3_E; num_steps = n , close_loop = false, face_type = :quad)
### ### ### FL_TOP 1st Layer, Merge All
VFL_1       =       [V7; V_N2; V_N3]
FFL_1       =       [F7;
                    [f .+ length(V7) for f in F_N2];
                    [f .+ length(V7) .+ length(V_N2) for f in F_N3]]
FFL_1, VFL_1=       mergevertices(FFL_1, VFL_1)
Eb          =       boundaryedges(FFL_1)     ## ??
indfix      =       unique(reduce(vcat, Eb)) ## ??
# nSmooth   =       5
# # λ       =       0.5                                   # Smoothening parameter
# # V           =       smoothmesh_laplacian(F, V, nSmooth, λ; constrained_points = indfix)
### ### ### FL_TOP FINAL EXTRUDE
EFL_T,VFL_T  =      extrudefaces(FFL_1, VFL_1; extent=tFL, direction=:positive, num_steps=3)
FFL_T        =      element2faces(EFL_T)
### ### ### WEB
VWeb_0       =      filter(p -> isapprox(p[3],0,atol=searchTol), VFL_T)
VWeb_S1      =      filter(p -> isapprox(p[1],-0.5*thcWeb,atol=searchTol), VWeb_0)
indSort      =      reverse(sortperm([v[2] for v in VWeb_S1]))
VWeb_S1      =      VWeb_S1[indSort]
VWeb_E1      =      [Point{3, Float64}(v[1], v[2], v[3]-depthWeb) for v in VWeb_S1]
FWeb1, VWeb1 =      loftlinear(VWeb_E1, VWeb_S1; num_steps = 10, close_loop = false, face_type = :quad)
EWeb2, VWeb2 =      extrudefaces(FWeb1, VWeb1; extent=thcWeb, direction=:negative, num_steps=3)
FWeb2        =      element2faces(EWeb2)
### ### ### FL_TOP 2nd Layer
VFL_2        =       [Point{3, Float64}(v[1], v[2], v[3]-depthWeb) for v in VFL_1]
FFL_2        =       FFL_1

EFL_B,VFL_B  =       extrudefaces(FFL_2, VFL_2; extent=tFL, direction=:positive, num_steps=3)
FFL_B        =       element2faces(EFL_B)

### Visualization
markersize = 25
linewidth  = 1
fig = Figure(size = (1200, 800))
# ax1 = Axis(fig[1, 1], xlabel = "X", ylabel = "Y")
ax1 = AxisGeom(fig[1, 1])

# scatter!(ax1, VFL_T     , markersize = markersize *0.25   , color = :Red)
# scatter!(ax1, VWeb_0    , markersize = markersize *0.25   , color = :Yellow)
# scatter!(ax1, VWeb_S1   , markersize = markersize *0.50   , color = :Green)
# scatter!(ax1, VWeb_E1   , markersize = markersize *0.50   , color = :Purple)
# scatter!(ax1, VFL_2     , markersize = markersize *0.50   , color = :Blue)

poly!(ax1, GeometryBasics.Mesh(VWeb2, FWeb2), strokewidth = linewidth, color = :white, strokecolor = :black, shading = FastShading, transparency = false)
# normalplot(ax1,FWeb2,VWeb2)
poly!(ax1, GeometryBasics.Mesh(VFL_B, FFL_B), strokewidth = linewidth, color = :white, strokecolor = :black, shading = FastShading, transparency = false)
# normalplot(ax1,FFL_B,VFL_B)
poly!(ax1, GeometryBasics.Mesh(VFL_T, FFL_T), strokewidth = linewidth, color = :white, strokecolor = :black, shading = FastShading, transparency = false)
# normalplot(ax1,FFL_T,VFL_T)

# wireframe!(ax1, GeometryBasics.Mesh(Vh, Eh), linewidth = 3, color = :red)
fig