using ISC
using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Statistics
using Comodo.Rotations
using Comodo.LinearAlgebra
using Printf

## INPUTS
r           =       1.5                     # Corner Radii
TH          =       15.0                    # Tooth Height      # at max    18      min ----, it shows error
TD          =       11.5                    # Tooth Depth       # at max    20      min ----, it shows error
TG          =       11.0                    # Tooth Gap         # at max    unl     min 11  ,
bFL         =       50.0                    # Width of FL       # at max    unl     min 34  , it shows error (technically no concern with optimization)
tFL         =       6.0                     # FL Thickness
beFL        =       bFL - 2 * TD            # Effective W.of FL
FL_B_L      =       30.0
nt_FL       =       2                       # because it was duplicating the value in 7-SP, which is included in this
nSections   =       nt_FL
disp        =       3.5                     # disaplacement applied U2 as BC
g           =       1.0                     # FL-FL Gap
searchTol   =       1e-6 
## Mesh Spacing Calculator
SP          =       1.5
num_cor_h   =       1                       # multiplier
smoothness  =       5

n           =       ceil(Int,((0.5*TH+r)/SP))
n_num1      =       n
n_num2      =       ceil(Int,(num_cor_h*((0.5*TH-r)/SP)))
n_num3      =       max(3, ceil(Int, ((0.5 * TG - r) / SP))) # it will make max. value of 3 or if formula gives more)
n_num4      =       floor(Int,((TD-r)/SP))
n_numN22    =       n_num3
n_numN3     =       ceil(Int,((FL_B_L-0.5*TH)/SP))
nThc        =       ceil(Int, tFL / SP) + 1 
## Shift Distances
shiftDistance = (TH + TG)
shift1        = 0.5*TG+TH+0.5*TG-0.5*g + g +  0.5*TG+TH+0.5*TG-0.5*g  
                                            # shift assembly T FL from center of B FL (y-direction only)
shift2x       = 0.5*bFL-TD+r                # shift assembly center of B FL to center of SP (x-direction only)
shift2y       = (nSections+(nSections+1))*(0.5*TG+0.5*TH)                     # shift assembly center of B FL to center of SP (y-direction only)    
######################################      ## 1st - Lofting from HI to AB 
e1_h        =       0.5*TH+r
e1_w        =       0.5*beFL-0.5*TH                                 
V1_S        =       Vector{Point{3, Float64}}(undef, n)          
v_1         =       Point{3, Float64}(0.0, 0.0 , 0.0)
v_2         =       Point{3, Float64}(0.0, e1_h , 0.0)
for (i, w) in enumerate(range(0.0, 1.0, n))              # Linear interpolation
    V1_S[i] =       (1.0 - w) * v_1 + w * v_2          
end
V1_E        =       [Point{3, Float64}(v[1]-e1_w, v[2], v[3]) for v in V1_S]
F1, V1      =       loftlinear(V1_S, V1_E; num_steps = max(3,n_num1-3), close_loop = false, face_type = :quad)
######################################      ## 2nd - Lofting from AB BC to XY
e2_w        =       0.5*TH 
Vr          =       [Point{3, Float64}(r * cosd(t), r * sind(t), 0.0) for t in range(0.0, 90, n+(n-1))]
Vr          =       [Point{3, Float64}(v[1]-e1_w-e1_h, v[2], v[3]) for v in Vr]
Vb          =       Vector{Point{3, Float64}}(undef, n)  
Vb1_E       =       [Point{3, Float64}(v[1], v[2], v[3]) for v in V1_E]
Vb2_E       =       Vector{Point{3, Float64}}(undef, n)          
v_1         =       Point{3, Float64}(-e1_w     , e1_h , 0.0)
v_2         =       Point{3, Float64}(-e1_w-e1_h, e1_h , 0.0)
for (i, w) in enumerate(range(0.0, 1.0, n))              # Linear interpolation
    Vb2_E[i] = (1.0 - w) * v_1 + w * v_2                
end
Vb          =       [Vb1_E;Vb2_E]
Vb          =       unique(Vb, dims=1)
F2, V2      =       loftlinear(Vb, Vr; num_steps = n_num2, close_loop = false, face_type = :quad)
######################################      ## 3rd - Lofting from HX to GF
e3_h        =       0.5*TG-r   
V3_S1       =       filter(p -> isapprox(p[2], 0.0 ,atol=searchTol), V1)
V3_S2       =       filter(p -> isapprox(p[2], 0.0 ,atol=searchTol), V2)
V3_S        =       [V3_S1;V3_S2]
V3_S        =       unique(V3_S, dims=1)
V3_E        =       [Point{3, Float64}(v[1], v[2]-e3_h, v[3]) for v in V3_S]
F3, V3      =       loftlinear(V3_S, V3_E; num_steps = n_num3, close_loop = false, face_type = :quad)
######################################      ## 4th - Lofting from CY to DE
e4_w        =       TD-r
V4_S        =       filter(p -> isapprox(p[1], -e1_w-e1_h ,atol=searchTol), V2)
V4_E        =       [Point{3, Float64}(v[1]-e4_w, v[2], v[3]) for v in V4_S]
F4, V4      =       loftlinear(V4_E, V4_S; num_steps = n_num4, close_loop = false, face_type = :quad)
######################################      ## merging 1234
V1234       =       [V1; V2; V3; V4]
F1234       =       [F1;
                    [f .+ length(V1) for f in F2];
                    [f .+ length(V1) .+ length(V2) for f in F3];
                    [f .+ length(V1) .+ length(V2) .+ length(V3) for f in F4]]
F1234, V1234=       mergevertices(F1234, V1234)
Eb1234      =       boundaryedges(F1234)                  # Constrained Smoothing 
V1234     .+=       Point{3, Float64}(0.0, +e3_h, 0.0 )
#########################################    ## Copy - 1 (Copy of basic element)
V_N1        =       [V1234;
                    [Point{3, Float64}(-v[1],  v[2], v[3]) for v in V1234];
                    [Point{3, Float64}( v[1], -v[2], v[3]) for v in V1234];
                    [Point{3, Float64}(-v[1], -v[2], v[3]) for v in V1234]]    
F_N1        =       [F1234;
                    [reverse(f)  .+    length(V1234)      for f in F1234];
                    [reverse(f)  .+    length(V1234)*2    for f in F1234];
                    [f           .+    length(V1234)*3    for f in F1234]]
#########################################    ## Copy - to create nSections
nf          =       length(F_N1)
nv          =       length(V_N1)
F_N1c       =       Vector{QuadFace{Int64}}(undef,nf*nSections)
V_N1c       =       Vector{Point{3,Float64}}(undef,nv*nSections)
i_f         =       1
i_v         =       1
s           =       0
for q in 1:nSections
    F_N1c[i_f:i_f+nf-1]  = [f.+s for f in F_N1]
    V_N1c[i_v: i_v+nv-1] = [Point{3, Float64}(v[1], v[2] - (q-1)*shiftDistance, v[3]) for v in V_N1]
    global i_f += nf
    global i_v += nv
    global s   += nv
end
F_N1c, V_N1c=      mergevertices(F_N1c, V_N1c)
#########################################    ## Loft Top-Part 1
V_N21       =       [[Point{3, Float64}(v[1],  -v[2]  .+ shiftDistance-0.5*TG+r, v[3]) for v in V1];
                    [Point{3, Float64}(v[1],   -v[2]  .+ shiftDistance-0.5*TG+r, v[3]) for v in V2];
                    [Point{3, Float64}(v[1],   -v[2]  .+ shiftDistance-0.5*TG+r, v[3]) for v in V4]]

F_N21       =       [[reverse(f) for f in F1];
                    [reverse(f) .+ length(V1) for f in F2]; 
                    [reverse(f) .+ length(V1) .+ length(V2) for f in F4]]
F_N21, V_N21=       mergevertices(F_N21, V_N21)
Eb_N21      =       boundaryedges(F_N21)
   
#########################################    ## Loft Top-Part 2
y           =       0.5*TG - 0.5*g - r
V_N22_S     =       filter(p -> isapprox(p[2], 0.5*TG + TH + r ,atol=searchTol), V_N21)
indSort     =       reverse(sortperm([v[1] for v in V_N22_S]))
V_N22_S     =       V_N22_S[indSort]
V_N22_E     =       [Point{3, Float64}(v[1],  v[2]  .+ y, v[3]) for v in V_N22_S]
indSort     =       reverse(sortperm([v[1] for v in V_N22_S]))
V_N22_E     =       V_N22_E[indSort]
F_N22, V_N22=       loftlinear(V_N22_E, V_N22_S ; num_steps = n_numN22, close_loop = false, face_type = :quad)
#########################################    ## Merging 21 & 22
V_N2122     =       [V_N21; V_N22]                    
F_N2122     =       [F_N21;
                    [f .+ length(V_N21) for f in F_N22]]
F_N2122, V_N2122 =  mergevertices(F_N2122, V_N2122)
Eb_N2122    =       boundaryedges(F_N2122)                        
#########################################    ## Merging Top Part 2
V_N2        =       [V_N2122;
                    [Point{3, Float64}(-v[1],  v[2], v[3]) for v in V_N2122]]     # top right    # bottom left
F_N2        =       [F_N2122;
                    [reverse(f)  .+    length(V_N2122)     for f in F_N2122]]
F_N2, V_N2  =       mergevertices(F_N2, V_N2)
Eb_N2       =       boundaryedges(F_N2) 
#########################################    ## Loft Bottom
V_N3_S      =       filter(p -> isapprox(p[2],-(nSections-0.5)*(TG+TH),atol=searchTol), V_N1c)
indSort     =       reverse(sortperm([v[1] for v in V_N3_S]))
V_N3_S      =       V_N3_S[indSort]
V_N3_E      =       [Point{3, Float64}(v[1], -(FL_B_L-0.5*TH)+v[2] , v[3]) for v in V_N3_S]
F_N3, V_N3  =       loftlinear(V_N3_S, V_N3_E; num_steps = n_numN3 , close_loop = false, face_type = :quad)

#########################################    ## Merge All
V           =       [V_N1c; V_N2; V_N3]
F           =       [F_N1c;
                    [f .+ length(V_N1c) for f in F_N2];
                    [f .+ length(V_N1c) .+ length(V_N2) for f in F_N3]]
F, V        =       mergevertices(F, V)
Eb          =       boundaryedges(F)  
indfix      =       unique(reduce(vcat, Eb))
nSmooth     =       smoothness
λ           =       0.5                                   # Smoothening parameter
V           =       smoothmesh_laplacian(F, V, nSmooth, λ; constrained_points = indfix)
#########################################    ## Extrude
Eh,Vh       =       extrudefaces(F, V; extent=tFL, direction=:positive, num_steps=nThc)
Fh          =       element2faces(Eh)
Fhb         =       boundaryfaces(Fh)

#########################################    ## Making SETS for ABAQUS INPUT FILE 
#########################################    ## make Set-1 for the BC fix  @ B FL Bot End, pts from Vh having y val = V_N3_E 
indext_in_V_N3_E = 2                         # targeting 2nd (y) value in each node of V_N4_E
target_y    =       V_N3_E[indext_in_V_N3_E][2]
V_Set_1     =       filter(p -> isapprox(p[2], target_y, atol=searchTol), Vh)
node_num_V_Set_1 = findall(p -> isapprox(p[2], target_y, atol=searchTol), Vh)

##################################################    ## Creating Input File (basic final file is file 1)
file1       =       joinpath(iscdir(), "assets", "temp", "temp_FL.inp")
file2       =       joinpath(iscdir(), "assets", "temp", "temp_SP.inp")
file3       =       joinpath(iscdir(), "assets", "inp", "Material_YG.inp")

open(file1,"w") do io
    println(io, "*Heading")
    println(io, "*** Job name: merged_file Model name: merged_file")
    println(io, "** Generated by: Abaqus/CAE 2023")
    println(io, "*Preprint, echo=NO, model=NO, history=NO, contact=NO")
    println(io, "**")
    println(io, "** PARTS")
    println(io, "**")
    println(io, "*Part, name=FL")
    println(io, "* Node")
    for (i,v) in enumerate(Vh)
        println(io,@sprintf("%i", i) * ", " * join([@sprintf("%.16e",x) for x in v],", "))
    end
    println(io,"*Element, type=C3D8R")
    for (i,e) in enumerate(Eh)
        println(io,@sprintf("%i", i) * ", " * join([@sprintf("%i",x) for x in e],", "))
    end
    println(io, "*Elset, elset=\"Set-5 (Full FL)\", generate")
    println(io, "1,  ", length(Eh), ",    1") # need-make 740 a function
    println(io, "** Section: ISC")
    println(io, "*Solid Section, elset=\"Set-5 (Full FL)\", material=YG-STEEL")
    println(io, ",")
    println(io, "*End Part") 
    println(io, "**")
    open(file2, "r") do f2 
        for line in eachline(f2)
            println(io, line)  # Write each line from the second file to the output
        end
        println(io, "*Elset, elset=\"Set-6 (Full SP)\", generate")
        println(io, "1,  ", length(EN6), ",    1")   ## need-make 740 a function
        println(io, "** Section: ISC")
        println(io, "*Solid Section, elset=\"Set-6 (Full SP)\", material=YG-STEEL")
        println(io, ",")

        println(io, "*End Part") 
        println(io, "**")
        println(io, "**")
    end
    println(io, "** ASSEMBLY")                       ## Assembly code for FL            
    println(io, "**")
    println(io, "*Assembly, name=Assembly")
    println(io, "**")
    println(io, "*Instance, name=B_FL, part=FL")                 
    println(io, "*End Instance")
    println(io, "**")
    println(io, "*Instance, name=T_FL, part=FL")
                                                    # distance from origin B_FL to inverted origin T_FL
    println(io, "           0.,            ", shift1,",           0.")     
    println(io, "           0.,            ", shift1,",           0.,            0.,            ", shift1,",           1.,         180.") 
    println(io, "*End Instance")                    
    println(io, "**")
    println(io, "*Instance, name=R_SP, part=SP")    ## Assembly code for SP              
    println(io, "           ", shift2x,",    ", shift2y,",            3.")     
    println(io, "           ", shift2x,",    ", shift2y,",            3.,           ", shift2x,",                ", -shift2y,",           3.,         -90.")
    println(io, "*End Instance")
    println(io, "**")
    println(io, "*Instance, name=L_SP, part=SP")
    println(io, "           ", -shift2x,",    ", shift2y,",           3.")     
    println(io, "           ", -shift2x,",    ", shift2y,",           3.,            ", -shift2x,",              ", -shift2y,",           3.,         90.")
    println(io, "*End Instance")
    println(io, "**")
    #### Creation of Node Set 1 (BC 1-Fixed Support at B FL) 
    println(io, "*Nset, nset=Set-1, instance=B_FL")
    chunk_size = 8
    for i in 1:chunk_size:length(node_num_V_Set_1)
    chunk = node_num_V_Set_1[i:min(i+chunk_size-1, end)]
    println(io, join(chunk, ", "))
    end
    #### Creation of Node Set 2 (BC 1-Disp. applied at T FL) 
    println(io, "*Nset, nset=Set-2, instance=T_FL")
    chunk_size = 8
    for i in 1:chunk_size:length(node_num_V_Set_1)
    chunk = node_num_V_Set_1[i:min(i+chunk_size-1, end)]
    println(io, join(chunk, ", "))
    end
    #### Creation of Node Set 3 (BC 1-Disp. applied at R_SP) 
    println(io, "*Nset, nset=Set-3, instance=R_SP")
    println(io, join(node_num_V_Set_3, ", "))
    #### Creation of Node Set 4 (BC 1-Disp. applied at L_SP) 
    println(io, "*Nset, nset=Set-4, instance=L_SP")
    println(io, join(node_num_V_Set_3, ", "))

    println(io, "*End Assembly")
    println(io, "**")
end
#### Defining Material (Seperate file of YG Model is pasted in the folder)
open(file1, "a") do io   
    open(file3, "r") do f3
        for line in eachline(f3)
            println(io, line)
        end
    end
    println(io, "** INTERACTION PROPERTIES")
    println(io, "**") 
    println(io, "*Surface Interaction, name=fric")
    println(io, "1.,")
    println(io, "*Friction, slip tolerance=0.005")
    println(io, "0.05,")
    println(io, "**")

    #### Boundary Condition
    println(io, "** BOUNDARY CONDITIONS")
    println(io, "**")
    println(io, "** Name: BC-1 Type: Displacement/Rotation")
    println(io, "*Boundary")
    println(io, "Set-1, 1, 1") 
    println(io, "** Name: BC-2 Type: Displacement/Rotation")
    println(io, "*Boundary")
    println(io, "Set-1, 2, 2")
    println(io, "** Name: BC-3 Type: Displacement/Rotation")
    println(io, "*Boundary")
    println(io, "Set-1, 3, 3")

    println(io, "**") 
    println(io, "** INTERACTIONS")
    println(io, "**") 
    println(io, "** Interaction: ALL")
    println(io, "*Contact")
    println(io, "*Contact Inclusions, ALL EXTERIOR")
    println(io, "*Contact Property Assignment")
    println(io, ",  , fric")

    println(io, "** ----------------------------------------------------------------")
    #### Step Load
    println(io, "**")
    println(io, "** STEP: Load")
    println(io, "**")
    println(io, "*Step, name=Load, nlgeom=YES")
    println(io, "Apply Load")
    println(io, "*Static, stabilize, factor=0.0002, allsdtol=0, continue=NO")
    println(io, "1., 1., 1e-05, 1.")
    # BC4
    println(io, "**")
    println(io, "** BOUNDARY CONDITIONS")
    println(io, "**")
    println(io, "** Name: BC-4 Type: Displacement/Rotation")
    println(io, "*Boundary")
    println(io, "Set-2, 1, 1") 
    println(io, "** Name: BC-5 Type: Displacement/Rotation")
    println(io, "*Boundary")
    println(io, "Set-2, 2, 2, ", disp) #later, you can apply disp. as a function
    println(io, "** Name: BC-6 Type: Displacement/Rotation")
    println(io, "*Boundary")
    println(io, "Set-2, 3, 3")
    println(io, "** Name: BC-7 Type: Displacement/Rotation")

    println(io, "*Boundary")
    println(io, "Set-3, 1, 1")

    println(io, "** Name: BC-8 Type: Displacement/Rotation")
    println(io, "*Boundary")
    println(io, "Set-3, 3, 3")

    println(io, "** Name: BC-9 Type: Displacement/Rotation")
    println(io, "*Boundary")
    println(io, "Set-4, 1, 1")

    println(io, "** Name: BC-10 Type: Displacement/Rotation")
    println(io, "*Boundary")
    println(io, "Set-4, 3, 3")

    println(io, "**")
    println(io,"** OUTPUT REQUESTS")
    println(io,"**")
    println(io,"*Restart, write, frequency=0")
    println(io,"**") 
    println(io,"** FIELD OUTPUT: F-Output-1")
    println(io,"**") 
    println(io,"*Output, field, variable=PRESELECT")
    println(io,"**") 
    println(io,"** HISTORY OUTPUT: H-Output-1")
    println(io,"**") 
    println(io,"*Output, history, variable=PRESELECT")
    println(io,"*End Step")
end
include("6_SP.jl")

#########################################    ## Visualization

Fhp, Vhp = separate_vertices(Fh, Vh)

markersize = 25
linewidth  = 1
fig = Figure(size = (1200, 800))
ax1 = AxisGeom(fig[1, 1], title="Flange design and mesh")
meshplot!(ax1, Fhp, Vhp, strokewidth = linewidth, color = :white, strokecolor = :black)
# normalplot(ax1,Fp,Vp)
fig