"""
    iscdir()

# Description 

This function simply returns the string for the ISC package path. This is helpful for
instance to load items, such as meshes, from the `assets` folder. 
"""
function iscdir()
    pkgdir(@__MODULE__)
end

"""
    mesh_sideplate(; B_SP=20.0, TG=11.0, TH=15.0, g = 0.5, TThc=6.0, tSP=6.0, mesh_spacing=4.5, nSections=2, searchTol = 1e-6)    

Creates mesh for sideplate

# Description 
fdasfasfjasfasfa f 
"""
function mesh_sideplate(; B_SP=20.0, TG=11.0, TH=15.0, g = 0.5, TThc=6.0, tSP=6.0, mesh_spacing=4.5, nSections=2, searchTol = 1e-6)    
    # Hole_H  = 2.0*g + TH
    Hole_W  = 2.0*g + TThc 
    nSectionsHoles = 2*nSections
    n = ceil(Int, B_SP / mesh_spacing)+1

    #################################################        ## 1st Lofting 
    e1_w    =   0.5*(B_SP-Hole_W)                            # 1st Element-Width

    V1_S    =   Vector{Point{3, Float64}}(undef, n)          # Create Points - r segment AB and BC
    v_1     =   Point{3, Float64}(e1_w, 0.0 , 0.0)
    v_2     =   Point{3, Float64}(0.0 , 0.0 , 0.0)
    for (i, w) in enumerate(range(0.0, 1.0, n))              # Linear interpolation
        V1_S[i] = (1.0 - w) * v_1 + w * v_2                  # from v_1 to v_2 A to B
    end
    V1_E    =   [Point{3, Float64}(v[1], v[2]+ 0.5* TH + g, v[3]) for v in V1_S]
    F1, V1  =   loftlinear(V1_E, V1_S; num_steps = n+1, close_loop = false, face_type = :quad)
    ###################################################      ## 2nd Lofting 
    V2_S    =   [Point{3, Float64}(v[1], v[2], v[3]) for v in V1_E]
    V2_E    =   [Point{3, Float64}(v[1], v[2] + 0.5*TG - g , v[3]) for v in V2_S]
    F2, V2  =   loftlinear(V2_E, V2_S; num_steps = n-1 , close_loop = false, face_type = :quad)
    ###################################################      ## 3rd Lofting 
    V3_S    =   filter(p -> isapprox(p[1],e1_w,atol=searchTol), V2) 
    indSort =   reverse(sortperm([v[1] for v in V3_S]))
    V3_E    =   [Point{3, Float64}(v[1]+ Hole_W, v[2], v[3]) for v in V3_S]
    F3, V3  =   loftlinear(V3_S, V3_E; num_steps = n, close_loop = false, face_type = :quad)
    ###################################################      ## 4th Lofting 
    V4_S    =   filter(p -> isapprox(p[1],e1_w + Hole_W,atol=searchTol), V3) 
    indSort =   reverse(sortperm([v[1] for v in V3]))
    V4_E    =   [Point{3, Float64}(v[1] + e1_w, v[2] , v[3]) for v in V4_S]
    F4, V4  =   loftlinear(V4_S, V4_E; num_steps = n, close_loop = false, face_type = :quad)
    ###################################################      ## 5th Lofting 
    V5_S    =   [Point{3, Float64}(v[1] + e1_w + Hole_W , v[2], v[3]) for v in V1_E]
    indSort =   reverse(sortperm([v[1] for v in V5_S]))
    V5_E    =   [Point{3, Float64}(v[1] + e1_w + Hole_W , v[2], v[3]) for v in V1_S]
    F5, V5  =   loftlinear(V5_S, V5_E; num_steps = n+1, close_loop = false, face_type = :quad)
    ################################################    ## Merge - EL 1 2 3 4 5
    V12345  = [V1; V2; V3; V4; V5]
    F12345  = [F1;
            [f .+ length(V1) for f in F2];
            [f .+ length(V1) .+ length(V2) for f in F3];
            [f .+ length(V1) .+ length(V2) .+ length(V3) for f in F4]
            [f .+ length(V1) .+ length(V2) .+ length(V3) .+ length(V4) for f in F5]]

    # ################################################    ## Copy - 1 (Copy of basic element)
    V12345  .+= Point{3, Float64}(-0.5*B_SP, 0.0 , 0.0 )

    VN1    = [V12345; 
            [Point{3, Float64}(v[1], -v[2], v[3]) for v in V12345]]   
    FN1    = [F12345; 
            [reverse(f)  .+    length(V12345)      for f in F12345]]
    # # ################################################    ## Copy - n 
    shiftDistance = TH + TG

    nv  = length(VN1)
    nf  = length(FN1)

    FN1c = Vector{QuadFace{Int64}}(undef,nf*nSectionsHoles)
    VN1c = Vector{Point{3,Float64}}(undef,nv*nSectionsHoles)
    i_f = 1
    i_v = 1
    s = 0
    for q in 1:nSectionsHoles
        FN1c[i_f:i_f+nf-1]  = [f.+s for f in FN1] 
        VN1c[i_v: i_v+nv-1] = [Point{3, Float64}(v[1], v[2] - (q-1)*shiftDistance, v[3]) for v in VN1] 
        i_f += nf 
        i_v += nv 
        s   += nv
    end
    FN1c, VN1c= mergevertices(FN1c, VN1c)
    ##################################################    ## Copy - 3 (Top)
    VN3_S     = filter(p -> isapprox(p[2],+0.5*(TG+TH),atol=searchTol), VN1c) 
    indSort = reverse(sortperm([v[1] for v in VN3_S]))
    VN3_S = VN3_S[indSort]
    VN3_E     = [Point{3, Float64}(v[1], v[2] + 0.5*TG - g, v[3]) for v in VN3_S]
    FN3, VN3  = loftlinear(VN3_E, VN3_S; num_steps = n, close_loop = false, face_type = :quad)
    ##################################################    ## Copy - 4 (Bottomp)
    VN4_S     = filter(p -> isapprox(p[2],-(nSectionsHoles-0.5)*(TG+TH),atol=searchTol), VN1c) 
    indSort = reverse(sortperm([v[1] for v in VN4_S]))
    VN4_S = VN4_S[indSort]
    VN4_E     = [Point{3, Float64}(v[1], v[2] - 0.5*TG + g, v[3]) for v in VN4_S]
    FN4, VN4  = loftlinear(VN4_S, VN4_E; num_steps = n, close_loop = false, face_type = :quad)
    ##################################################    ## Final (2D) - 5 (n + 3 + 4 - Copy & Merge)
    VN5 = [VN1c; VN3; VN4]
    FN5 = [FN1c;
            [f .+ length(VN1c) for f in FN3];
            [f .+ length(VN1c) .+ length(VN3) for f in FN4]]

    FN5, VN5 = mergevertices(FN5, VN5)
    ##################################################    ## Final (3D) Extrude
    EN6,VN6 = extrudefaces(FN5, VN5; extent=tSP, direction=:positive, num_steps=n)

    return EN6, VN6
end

function mesh_flange(; r =1.0, g =1.0 ,Gap = 1.0  ,TH = 10.0   ,TD =  11.5,TG =  11.0,bFL = 50.0, tFL = 6.0, nt_FL =  3 ,L =  98.0 ,searchTol=  1e-6)  
    ### DERIVED INPUTS
    shiftDist   =       (TH + TG)
    FL_B_L      =       L-(nt_FL *(shiftDist) + 0.5*TG - 0.5*Gap)
    beFL        =       bFL - 2 * TD            # Effective W.of FL
    nSections   =       nt_FL
    n = 5                                       # n=ceil(Int,((0.5*TH+r)/1.5))

    ### 1st LOFTING - EL 1 - HI to AB
    e1_h        =       0.5*TH+r
    e1_w        =       0.5*beFL-0.5*TH                                 
    V1_S        =       Vector{Point{3, Float64}}(undef, n)          
    v_1         =       Point{3, Float64}(0.0, 0.0 , 0.0)
    v_2         =       Point{3, Float64}(0.0, e1_h , 0.0)
    for (i, w) in enumerate(range(0.0, 1.0, n))             
        V1_S[i] =       (1.0 - w) * v_1 + w * v_2          
    end
    V1_E        =       [Point{3, Float64}(v[1]-e1_w, v[2], v[3]) for v in V1_S]
    F1, V1      =       loftlinear(V1_S, V1_E; num_steps = 10, close_loop = false, face_type = :quad)
    ### 2nd LOFTING - EL 2 - AB BC to XY       
    e2_w        =       0.5*TH 
    Vr          =       [Point{3, Float64}(r * cosd(t), r * sind(t), 0.0) for t in range(0.0, 90, n+(n-1))]
    Vr          =       [Point{3, Float64}(v[1]-e1_w-e1_h, v[2], v[3]) for v in Vr]
    Vb          =       Vector{Point{3, Float64}}(undef, n)  
    Vb1_E       =       [Point{3, Float64}(v[1], v[2], v[3]) for v in V1_E]
    Vb2_E       =       Vector{Point{3, Float64}}(undef, n)          
    v_1         =       Point{3, Float64}(-e1_w     , e1_h , 0.0)
    v_2         =       Point{3, Float64}(-e1_w-e1_h, e1_h , 0.0)
    for (i, w) in enumerate(range(0.0, 1.0, n))              
        Vb2_E[i] = (1.0 - w) * v_1 + w * v_2                
    end
    Vb          =       [Vb1_E;Vb2_E]
    Vb          =       unique(Vb, dims=1)
    F2, V2      =       loftlinear(Vb, Vr; num_steps = 10, close_loop = false, face_type = :quad)
    ### 3rd LOFTING - EL 3 - HX to AF
    e3_h        =       0.5*TG-r   
    V3_S1       =       filter(p -> isapprox(p[2], 0.0 ,atol=searchTol), V1)
    V3_S2       =       filter(p -> isapprox(p[2], 0.0 ,atol=searchTol), V2)
    V3_S        =       [V3_S1;V3_S2]
    V3_S        =       unique(V3_S, dims=1)
    V3_E        =       [Point{3, Float64}(v[1], v[2]-e3_h, v[3]) for v in V3_S]
    F3, V3      =       loftlinear(V3_S, V3_E; num_steps = 5, close_loop = false, face_type = :quad)
    ### 4th LOFTING - EL 4 - YC to DE
    e4_w        =       TD-r
    V4_S        =       filter(p -> isapprox(p[1], -e1_w-e1_h ,atol=searchTol), V2)
    V4_E        =       [Point{3, Float64}(v[1]-e4_w, v[2], v[3]) for v in V4_S]
    F4, V4      =       loftlinear(V4_E, V4_S; num_steps = n, close_loop = false, face_type = :quad)
    ### Merging 1234
    V1234       =       [V1; V2; V3; V4]
    F1234       =       [F1;
                        [f .+ length(V1) for f in F2];
                        [f .+ length(V1) .+ length(V2) for f in F3];
                        [f .+ length(V1) .+ length(V2) .+ length(V3) for f in F4]]
    F1234, V1234=       mergevertices(F1234, V1234)
    Eb1234      =       boundaryedges(F1234)                  
    V1234     .+=       Point{3, Float64}(0.0, +e3_h, 0.0)
    ### Copy - 1 (Copy of basic element)
    V_N1        =       [V1234;
                        [Point{3, Float64}(-v[1],  v[2], v[3]) for v in V1234];
                        [Point{3, Float64}( v[1], -v[2], v[3]) for v in V1234];
                        [Point{3, Float64}(-v[1], -v[2], v[3]) for v in V1234]]    
    F_N1        =       [F1234;
                        [reverse(f)  .+    length(V1234)      for f in F1234];
                        [reverse(f)  .+    length(V1234)*2    for f in F1234];
                        [f           .+    length(V1234)*3    for f in F1234]]
    ### Copy - nSections
    nf          =       length(F_N1)
    nv          =       length(V_N1)
    F_N1c       =       Vector{QuadFace{Int64}}(undef,nf*nSections)
    V_N1c       =       Vector{Point{3,Float64}}(undef,nv*nSections)
    i_f         =       1
    i_v         =       1
    s           =       0
    for q in 1:nSections
        F_N1c[i_f:i_f+nf-1]  = [f.+s for f in F_N1]
        V_N1c[i_v: i_v+nv-1] = [Point{3, Float64}(v[1], v[2] - (q-1)*shiftDist, v[3]) for v in V_N1]
        i_f += nf
        i_v += nv
        s   += nv
    end
    F_N1c, V_N1c=      mergevertices(F_N1c, V_N1c)
    ### Loft Top-Part 21 (Left)
    V_N21       =       [[Point{3, Float64}(v[1],  -v[2]  .+ shiftDist-0.5*TG+r, v[3]) for v in V1];
                        [Point{3, Float64}(v[1],   -v[2]  .+ shiftDist-0.5*TG+r, v[3]) for v in V2];
                        [Point{3, Float64}(v[1],   -v[2]  .+ shiftDist-0.5*TG+r, v[3]) for v in V4]]

    F_N21       =       [[reverse(f) for f in F1];
                        [reverse(f) .+ length(V1) for f in F2]; 
                        [reverse(f) .+ length(V1) .+ length(V2) for f in F4]]
    F_N21, V_N21=       mergevertices(F_N21, V_N21)
    ### Loft Top-Part 22 (Left)
    y           =       0.5*TG - 0.5*g - r
    V_N22_S     =       filter(p -> isapprox(p[2], 0.5*TG + TH + r ,atol=searchTol), V_N21)
    indSort     =       reverse(sortperm([v[1] for v in V_N22_S]))
    V_N22_S     =       V_N22_S[indSort]
    V_N22_E     =       [Point{3, Float64}(v[1],  v[2]  .+ y, v[3]) for v in V_N22_S]
    indSort     =       reverse(sortperm([v[1] for v in V_N22_S]))
    V_N22_E     =       V_N22_E[indSort]
    F_N22, V_N22=       loftlinear(V_N22_E, V_N22_S ; num_steps = 10, close_loop = false, face_type = :quad)
    ### Loft Top-Part Merge 21 & 22
    V_N2122     =       [V_N21; V_N22]                    
    F_N2122     =       [F_N21;
                        [f .+ length(V_N21) for f in F_N22]]
    F_N2122, V_N2122 =  mergevertices(F_N2122, V_N2122)
    Eb_N2122    =       boundaryedges(F_N2122)        #??
    ### Loft Top-Part Merge & Mirror 
    V_N2        =       [V_N2122;
                        [Point{3, Float64}(-v[1],  v[2], v[3]) for v in V_N2122]]     # top right    # bottom left
    F_N2        =       [F_N2122;
                        [reverse(f)  .+    length(V_N2122)     for f in F_N2122]]
    F_N2, V_N2  =       mergevertices(F_N2, V_N2)
    Eb_N2       =       boundaryedges(F_N2)
    ### Loft Bottom-Part
    V_N3_S      =       filter(p -> isapprox(p[2],-(nSections-0.5)*(TG+TH),atol=searchTol), V_N1c)
    indSort     =       reverse(sortperm([v[1] for v in V_N3_S]))
    V_N3_S      =       V_N3_S[indSort]
    V_N3_E      =       [Point{3, Float64}(v[1], -(FL_B_L-0.5*TH)+v[2] , v[3]) for v in V_N3_S]
    F_N3, V_N3  =       loftlinear(V_N3_S, V_N3_E; num_steps = n , close_loop = false, face_type = :quad)
    ### Merge All
    V           =       [V_N1c; V_N2; V_N3]
    F           =       [F_N1c;
                        [f .+ length(V_N1c) for f in F_N2];
                        [f .+ length(V_N1c) .+ length(V_N2) for f in F_N3]]
    F, V        =       mergevertices(F, V)
    Eb          =       boundaryedges(F)  
    indfix      =       unique(reduce(vcat, Eb))
    nSmooth     =       5
    λ           =       0.5                                   # Smoothening parameter
    V           =       smoothmesh_laplacian(F, V, nSmooth, λ; constrained_points = indfix)
    ### Extrude All
    Eh,Vh       =       extrudefaces(F, V; extent=tFL, direction=:positive, num_steps=3)
 return Eh,Vh
end

function f1(x)
    return x + 1
end

function f2(x)
    return x + 2
end