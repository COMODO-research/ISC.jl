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