function comp_bolt(; 
                    Center              =           Point{3, Float64}(0.0 , 0.0 , 0.0),
                    r_BHead             =           6.0,
                    Th_BHead            =           5.3,
                    r_Thread            =           4.0,
                    L_Thread            =           12.4,
                    n                   =           3,)
        # =============================================================================
        # SCRIPT - START
        # =============================================================================
        # INPUTS - DERIVED
        
        n_Edge1             =           n

        # CREATE - CIRCLE 1 - Thread
        F_Cir1, V_Cir1      =           create_circle(r = r_Thread, Center = Center, n_Edge1 = n)
        Eb                  =           boundaryedges(F_Cir1)
        indBoundary         =           edges2curve(Eb; remove_last = true)
        V_Cir1_Bnd          =           V_Cir1[indBoundary]
        V_Cir1_Bnd          =           sort(V_Cir1_Bnd; by = p -> mod(atan(p[2] - Center[2], p[1] - Center[1]), 2π))

        # CREATE - CIRCLE 2 - BHead
        V_Cir2_Bnd          =           [Point{3, Float64}(r_BHead * sind(t), r_BHead * cosd(t), 0.0) for t in range(0.0, 360.0, length(V_Cir1_Bnd)+1)][1:end-1]
        V_Cir2_Bnd          =           sort(V_Cir2_Bnd; by = p -> mod(atan(p[2] - Center[2], p[1] - Center[1]), 2π))
        F_Cir2, V_Cir2      =           loftlinear(V_Cir2_Bnd, V_Cir1_Bnd; num_steps = n_Edge1, close_loop = true, face_type = :quad,)
        F_Cir2, V_Cir2      =           mergevertices(F_Cir2, V_Cir2)
        # dup_nodes         =          count_duplicate_nodes(V_Cir2 )
        # println("Total duplicate nodes in V_Cir2  = $dup_nodes")

        # CREATE - CIRCLE 3 - CIRCLE 1 + CIRCLE 2
        F_Cir3, V_Cir3      =           combine_faces_nodes([V_Cir1, V_Cir2], [F_Cir1, F_Cir2], [:sam, :sam],)
        F_Cir3, V_Cir3      =           mergevertices(F_Cir3, V_Cir3)

        # EXTRUDE 1 
        E_Out, B_Out        =           extrudefaces(F_Cir3, V_Cir3; extent=Th_BHead, direction=:negative, num_steps=4)
        A_Out               =           element2faces(E_Out)

        # EXTRUDE 2 
        E_In, B_In          =           extrudefaces(F_Cir1, V_Cir1; extent=(0.5*L_Thread), direction=:positive, num_steps=5)
        A_In                =           element2faces(E_In)

        # BOLT _ HALF 1
        B_BHalf1            =           [B_Out; B_In]
        E_BHalf1            =           [E_Out;
                                        [e .+ length(B_Out) for e in E_In]]
        B_BHalf1, E_BHalf1  =           merge_nodes_and_update_elements(B_BHalf1, E_BHalf1)
        A_BHalf1            =           element2faces(E_BHalf1)
        B_BHalf1            =           [Point{3, Float64}(v[1], v[2], v[3]-(0.5*L_Thread)) for v in B_BHalf1]

        # BOLT - HALF 2
        B_BHalf2            =           [Point{3, Float64}(v[1], v[2], 0.0-v[3]) for v in B_BHalf1]
        E_BHalf2            =           [typeof(e)(e[5], e[6], e[7], e[8], e[1], e[2], e[3], e[4]) for e in E_BHalf1]   
        A_BHalf2            =           element2faces(E_BHalf2)

        # BOLT - FULL - HALF 1 + HALF 2 
        B_Bolt              =           [B_BHalf1; B_BHalf2]
        E_Bolt              =           [E_BHalf1;
                                        [(e) .+ length(B_BHalf1) for e in E_BHalf2]]
        B_Bolt, E_Bolt      =           merge_nodes_and_update_elements(B_Bolt, E_Bolt)
        A_Bolt              =           element2faces(E_Bolt)
        A_Bolt_Bnd          =           boundaryfaces(E_Bolt)

        # BOLT - FULL - ROTATE - CW 90 DEG 
        B_Bolt              =           [Point{3, Float64}(v[3] , v[2], -v[1]) for v in B_Bolt]
        # =============================================================================
        # SCRIPT - END
        # =============================================================================
    return A_Bolt_Bnd, B_Bolt, E_Bolt
end

function comp_WebStiffener(;
                    Wd_WS       = 20.0,
                    Th_WS       = 3.0 ,
                    Ht_WS       = 40.0, 
                    Mesh_Sp     = 20,
                    n_Th        = 3,)

        # =============================================================================
        # SCRIPT - START
        # =============================================================================
        # INPUTS - DERIVED
        nB          = ceil(Int, Wd_WS / Mesh_Sp)
        nL          = ceil(Int, Th_WS / Mesh_Sp)
        nH          = ceil(Int, Ht_WS / Mesh_Sp)

        n = nB+1
        # LOFT 1 -SURFACE
        e1_w        =   Wd_WS                          
        V1_S        =   Vector{Point{3, Float64}}(undef, n)        
        v_1         =   Point{3, Float64}(e1_w, 0.0 , 0.0)
        v_2         =   Point{3, Float64}(0.0 , 0.0 , 0.0)
        for (i, w) in enumerate(range(0.0, 1.0, n))              
            V1_S[i] = (1.0 - w) * v_1 + w * v_2                 
        end
        V1_E        =   [Point{3, Float64}(v[1], v[2] + Th_WS, v[3]) for v in V1_S]
        FWs, VWs    =   loftlinear(V1_E, V1_S; num_steps = n_Th, close_loop = false, face_type = :quad)

        # FINAL - EXTRUDE
        EWs,VWs     =   extrudefaces(FWs, VWs; extent= Ht_WS, direction=:positive, num_steps = nH)
        FWs         =   element2faces(EWs)
        VWs         =   [Point{3, Float64}(v[1] - 0.5*Wd_WS, v[2] - 0.5*Th_WS, v[3] - 0.5*Ht_WS) for v in VWs]

        ### CREATE SET - Web Stiffener for FL - R. Side
        VWs_SET_1   = [p for p in VWs
            if  p[1] >= +0.5*Wd_WS ||
            p[3] >= +0.5*Ht_WS  ||
            p[3] <= -0.5*Ht_WS ]

        VWs_SET_1_idx = findall(
            p -> p[1] >= +0.5*Wd_WS ||
                p[3] >= +0.5*Ht_WS ||
                p[3] <= -0.5*Ht_WS, VWs)

        ### CREATE SET - Web Stiffener for FL - L. Side
        VWs_SET_2 = [
            p for p in VWs
            if p[1] <= -0.5*Wd_WS ||
            p[3] >= +0.5*Ht_WS  ||
            p[3] <= -0.5*Ht_WS ]

        VWs_SET_2_idx = findall(
            p -> p[1] <= -0.5*Wd_WS ||
                p[3] >= +0.5*Ht_WS ||
                p[3] <= -0.5*Ht_WS, VWs)

        # =============================================================================
        # SCRIPT - START
        # =============================================================================
    return EWs, VWs, FWs, VWs_SET_1, VWs_SET_1_idx, VWs_SET_2, VWs_SET_2_idx 
end

function comp_SidePlate(; 
        Th_FL               =   6.0,                  
        nt_FL               =   2,                   
        TW                  =   10.0,                 
        TG                  =   11.0,                    
        g                   =   1.0,                 
        mesh_sp             =   1.0,
        b_SP                =   20.0,                
        Th_SP               =   6.0,
        searchTol = 1e-6)                  

        # =============================================================================
        # SCRIPT - START
        # =============================================================================
        # DERIVED INPUTS
        nSecHoles           =       2*nt_FL              
        Hole_B              =       0.5*g + Th_FL + 0.5*g    
        Hole_L              =       0.5*g + TW    + 0.5*g
        L1                  =       (0.5*TG -0.5*g)
        L2                  =       (0.5*TG + TW + 0.5*TG)     
        L_SP                =       L1 + nSecHoles*L2 + L1
        n_wd                =       ceil(Int, L_SP / mesh_sp)
        n_ht                =       ceil(Int, b_SP / mesh_sp)


        A1, B1              =       create_rect(;wd = L_SP, ht = b_SP, n_wd = n_wd, n_ht = n_ht, center = (0.0, 0.0, 0.0), plane = :YZ,)
        A2, B2              =       deepcopy((A1, B1))
        for i in 1:nSecHoles    
            y_value         =       (0.0 - 0.5*L_SP) + L1 + 0.5 * L2 + (i-1)*(TW + TG)
             A2, B2   =       cut_rect(A2, B2; Center = (0.0, y_value, 0.0), Rec1_w = Hole_L, Rec1_h = Hole_B, offset = 5.0, n_off = 4, plane = :YZ, Edge=:I,)
        end
        ESP, BSP            =       extrudefaces(A2, B2; extent=Th_SP, direction=:positive, num_steps=4)
        BSP                 =       [Point{3, Float64}(p[1] + 0.5* Th_SP, p[2], p[3]) for p in BSP]
        ASP                 =       element2faces(ESP)
        ASP_Bnd             =       boundaryfaces(ASP)

        # CREATE SET -  SP CORNER 
        BSP_SET_1           =       [p for p in BSP if
                                    (isapprox(p[1],  0.5*Th_SP, atol=searchTol) ||
                                    isapprox(p[1], -0.5*Th_SP, atol=searchTol)) &&
                                    (isapprox(p[2],  0.5*L_SP, atol=searchTol) ||
                                    isapprox(p[2], -0.5*L_SP, atol=searchTol)) &&
                                    (isapprox(p[3],  0.5*b_SP, atol=searchTol)  ||
                                    isapprox(p[3], -0.5*b_SP, atol=searchTol))]

        BSP_SET_1_idx       =       findall(p ->
                                    (isapprox(p[1],  0.5*Th_SP, atol=searchTol) ||
                                    isapprox(p[1], -0.5*Th_SP, atol=searchTol)) &&
                                    (isapprox(p[2],  0.5*L_SP, atol=searchTol) ||
                                    isapprox(p[2], -0.5*L_SP, atol=searchTol)) &&
                                    (isapprox(p[3],  0.5*b_SP, atol=searchTol)  ||
                                    isapprox(p[3], -0.5*b_SP, atol=searchTol)), BSP)
        # =============================================================================
        # SCRIPT - END
        # =============================================================================
    return ASP_Bnd, BSP, ESP, BSP_SET_1, BSP_SET_1_idx 
end

function comp_FL_v1(;
            TotalLength         = 197.0,                # Total length               - Beam to Beam
            Gap_BtoB            = 1.0,                  # Gap                        - Flange to Flange / Beam to Beam
            Edge_Support        = 50.0,                # Edge Dist. to support      - Beam to Beam
            Th_FL               = 6.0,                  # Thickness                  - Flange
            Th_Web              = 6.0,                  # Thickness                  - Web
            nt_FL               = 2,                    # no. of teeth (each side)   - Flange 
            b_FL                = 50.0,                 # Width                      - Flange   
            Ht_Sec              = 50.0,                 # Height                     - Section
            TW                  = 10.0,                 # Tooth Width                - Flange
            TG                  = 11.0,                 # Tooth Gap                  - Flange  
            TD                  = 11.5,                 # Tooth Depth                - Flange
            r                   = 1.0,                  # Corner Radii               - Flange   
            g                   = 1.0,                  # Tolerance                  - b/w Tooth and Holes
            searchTol           = 1e-6,)                 # Search Tolerance           - General

            # =============================================================================
            # SCRIPT - END
            # =============================================================================
            ### INPUTS
            Ht_Web          =       Ht_Sec-2*Th_FL         
            L               =       0.5*TotalLength-0.5*Gap_BtoB                    # Length of the Part excl. Gap_BtoB
            FL_B_L          =       L-(nt_FL*(TW + TG) + 0.5*TG - 0.5*Gap_BtoB)     # Bottom Length   - after the teeth section
            be_FL           =       b_FL - 2 * TD                                   # Effective Width - Flange
            nSections       =       nt_FL
            Pt_Load         =       TotalLength/6
            Pt_Support      =       L-Edge_Support
            searchTol       =       1e-6

            ### MESH CONTROL
            n_0             =       4                       
            n_1             =       4
            n_2             =       7                                                 # to link it
            n_3             =       5
            n_4             =       5
            n_5             =       n_0
            n_10            =       3
            n_11            =       n_1
            n_14            =       30
            n_16            =       4                                                # along Th_FL
            n_18            =       50  
            n_19            =       5   #   n_1+2                                            # along Th_Web

            ### FL_TOP - LOFT # 1, EL # 1 - HI to ED           
            e0_h            =       0.5*TW+r
            e0_w            =       0.5*Th_Web
            V0_S            =       Vector{Point{3, Float64}}(undef, n_0)
            v_1             =       Point{3, Float64}(0.0, 0.0 , 0.0)
            v_2             =       Point{3, Float64}(0.0, e0_h , 0.0)
            for (i, w) in enumerate(range(0.0, 1.0, n_0))             
                V0_S[i]     =       (1.0 - w) * v_1 + w * v_2          
            end
            V0_E            =       [Point{3, Float64}(v[1]-e0_w, v[2], v[3]) for v in V0_S]
            A1, B1          =       loftlinear(V0_S, V0_E; num_steps = n_1, close_loop = false, face_type = :quad)
            ### FL_TOP - LOFT # 2, EL # 2 - HI to AB
            e1_h            =       0.5*TW+r
            e1_w            =       0.5*be_FL-0.5*TW-e0_w  
            V1_S            =       [Point{3, Float64}(v[1], v[2], v[3]) for v in V0_E]                               
            V1_E            =       [Point{3, Float64}(v[1]-e1_w, v[2], v[3]) for v in V1_S]
            A2, B2          =       loftlinear(V1_S, V1_E; num_steps = n_2, close_loop = false, face_type = :quad)

            ### FL_TOP - LOFT # 3, EL # 3 - AB BC to XY       
            e2_w            =       0.5*TW 
            Vr              =       [Point{3, Float64}(r * cosd(t), r * sind(t), 0.0) for t in range(0.0, 90, n_0+(n_0-1))]
            Vr              =       [Point{3, Float64}(v[1]-e0_w-e1_w-e1_h, v[2], v[3]) for v in Vr]
            Vb              =       Vector{Point{3, Float64}}(undef, n_0)  
            Vb1_E           =       [Point{3, Float64}(v[1], v[2], v[3]) for v in V1_E]
            Vb2_E           =       Vector{Point{3, Float64}}(undef, n_0)          
            v_1             =       Point{3, Float64}(-e0_w-e1_w        , e1_h , 0.0)
            v_2             =       Point{3, Float64}(-e0_w-e1_w-e1_h   , e1_h , 0.0)
            for (i, w) in enumerate(range(0.0, 1.0, n_0))              
                Vb2_E[i]    = (1.0 - w) * v_1 + w * v_2                
            end
            Vb              =       [Vb1_E;Vb2_E]
            Vb              =       unique(Vb, dims=1)
            A3, B3          =       loftlinear(Vb, Vr; num_steps = n_3, close_loop = false, face_type = :quad)

            ### FL_TOP - LOFT # 4, EL # 4 - HX to GF
            e3_h            =       0.5*TG-r   
            V3_S0           =       filter(p -> isapprox(p[2], 0.0 ,atol=searchTol), B1)
            V3_S1           =       filter(p -> isapprox(p[2], 0.0 ,atol=searchTol), B2)
            V3_S2           =       filter(p -> isapprox(p[2], 0.0 ,atol=searchTol), B3)
            V3_S            =       [V3_S0;V3_S1;V3_S2]
            V3_S            =       unique(V3_S, dims=1)
            V3_E            =       [Point{3, Float64}(v[1], v[2]-e3_h, v[3]) for v in V3_S]
            A4, B4          =       loftlinear(V3_S, V3_E; num_steps = n_4, close_loop = false, face_type = :quad)

            ### FL_TOP - LOFT # 5, EL # 5 - YC to DE
            e4_w            =       TD-r
            V4_S            =       filter(p -> isapprox(p[1], -e0_w-e1_w-e1_h ,atol=searchTol), B3)
            V4_E            =       [Point{3, Float64}(v[1]-e4_w, v[2], v[3]) for v in V4_S]
            A5, B5          =       loftlinear(V4_E, V4_S; num_steps = n_5, close_loop = false, face_type = :quad)

            ### FL_TOP - LOFT # 6 (1,2,3,4,5)
            B6              =       [B1; B2; B3; B4; B5]
            A6              =       [A1;
                                    [f .+ length(B1) for f in A2];
                                    [f .+ length(B1) .+ length(B2) for f in A3];
                                    [f .+ length(B1) .+ length(B2) .+ length(B3) for f in A4];
                                    [f .+ length(B1) .+ length(B2) .+ length(B3).+ length(B4) for f in A5]]
            A6, B6          =       mergevertices(A6, B6)
            # Eb5           =       boundaryedges(A6)                  
            B6            .+=       Point{3, Float64}(0.0, +e3_h, 0.0) # if i turn it off, why code does not work

            ### FL_TOP - LOFT # 7 (BASIC ELEMENT)
            B7              =       [B6;
                                    [Point{3, Float64}(-v[1],  v[2], v[3]) for v in B6];
                                    [Point{3, Float64}( v[1], -v[2], v[3]) for v in B6];
                                    [Point{3, Float64}(-v[1], -v[2], v[3]) for v in B6]]    
            A7              =       [A6;
                                    [reverse(f)  .+    length(B6)      for f in A6];
                                    [reverse(f)  .+    length(B6)*2    for f in A6];
                                    [f           .+    length(B6)*3    for f in A6]]

            ### FL_TOP - LOFT # 8 (n COPY OF BASIC ELEMENTS)                        
            nf              =       length(A7)
            nv              =       length(B7)
            A8              =       Vector{QuadFace{Int64}}(undef,nf*nSections)
            B8              =       Vector{Point{3,Float64}}(undef,nv*nSections)
            i_f             =       1
            i_v             =       1
            s               =       0
            for q in 1:nSections
                A8[i_f:i_f+nf-1]  = [f.+s for f in A7]
                B8[i_v: i_v+nv-1] = [Point{3, Float64}(v[1], v[2] - (q-1)*(TW + TG), v[3]) for v in B7]
                    i_f += nf
                    i_v += nv
                    s   += nv
            end
            A8, B8          =       mergevertices(A8, B8)

            ### FL_TOP -  LOFT # 9 (TOP LEFT 1) - PMcG SCALED 
            B9              =       [[Point{3, Float64}(v[1],  -v[2]  .+ (TW + TG)-0.5*TG+r, v[3]) for v in B2];
                                    [Point{3, Float64}(v[1],   -v[2]  .+ (TW + TG)-0.5*TG+r, v[3]) for v in B3];
                                    [Point{3, Float64}(v[1],   -v[2]  .+ (TW + TG)-0.5*TG+r, v[3]) for v in B5]]

            A9              =       [[reverse(f) for f in A2];
                                    [reverse(f) .+ length(B2) for f in A3]; 
                                    [reverse(f) .+ length(B2) .+ length(B3) for f in A5]]
            A9, B9          =       mergevertices(A9, B9)

            ### FL_TOP -  LOFT # 10 (TOP LEFT 2) - PMcG SCALED 
            B10_S           =       filter(p -> isapprox(p[2], 0.5*TG + TW + r ,atol=searchTol), B9)
            indSort         =       reverse(sortperm([v[1] for v in B10_S]))
            B10_S           =       B10_S[indSort]
            B10_E           =       [Point{3, Float64}(v[1],  v[2]  .+ (0.5*TG - 0.5*Gap_BtoB - r), v[3]) for v in B10_S]
            indSort         =       reverse(sortperm([v[1] for v in B10_S]))
            B10_E           =       B10_E[indSort]
            A10,B10         =       loftlinear(B10_E, B10_S ; num_steps = n_10, close_loop = false, face_type = :quad)

            ### FL_TOP -  LOFT # 11 (TOP LEFT 3) - PMcG SCALED
            B11_S1          =       filter(p -> isapprox(p[1], -e0_w ,atol=searchTol), B9)
            B11_S2          =       filter(p -> isapprox(p[1], -e0_w ,atol=searchTol), B10)
            B11_S           =       [B11_S1;B11_S2]
            indSort         =       reverse(sortperm([v[2] for v in B11_S]))
            B11_S           =       B11_S[indSort]
            B11_S           =       unique(B11_S, dims=1)
            B11_E           =       [Point{3, Float64}(v[1]+e0_w, v[2], v[3]) for v in B11_S]
            indSort         =       reverse(sortperm([v[2] for v in B11_E]))
            B11_E           =       B11_E[indSort]
            A11,B11         =       loftlinear(B11_S, B11_E; num_steps = n_11, close_loop = false, face_type = :quad)

            ### FL_TOP -  LOFT # 12 (9 10 11) - PMcG SCALED
            B12             =       [B9; B10; B11]                    
            A12             =       [A9;
                                    [f .+ length(B9) for f in A10];
                                    [f .+ length(B9) .+ length(B10) for f in A11]]
            A12, B12        =       mergevertices(A12, B12)
            # Eb11          =       boundaryedges(A12)  

            ### FL_TOP -  LOFT # 13 (12 12) - PMcG SCALED
            B13             =       [B12;
                                    [Point{3, Float64}(-v[1],  v[2], v[3]) for v in B12]]     # top right    # bottom left
            A13             =       [A12;
                                    [reverse(f)  .+    length(B12)     for f in A12]]
            A13, B13        =       mergevertices(A13, B13)
            # Eb_N2         =       boundaryedges(A13)

            ### FL_TOP - LOFT # 14 - BOTTOM PART
            B14_S           =       filter(p -> isapprox(p[2],-(nSections-0.5)*(TG+TW),atol=searchTol), B8)
            indSort         =       reverse(sortperm([v[1] for v in B14_S]))
            B14_S           =       B14_S[indSort]
            B14_E           =       [Point{3, Float64}(v[1], -(FL_B_L-0.5*TW)+v[2] , v[3]) for v in B14_S]
            A14, B14        =       loftlinear(B14_S, B14_E; num_steps = n_14, close_loop = false, face_type = :quad)

            ### FL_TOP - LOFT # 15 
            B15             =       [B8; B13; B14]
            A15             =       [A8;
                                    [f .+ length(B8) for f in A13];
                                    [f .+ length(B8) .+ length(B13) for f in A14]]
            A15, B15        =       mergevertices(A15, B15)

            ### FL_TOP - LOFT # 16 - FINAL EXTRUDE
            E16,B16         =       extrudefaces(A15, B15; extent=Th_FL, direction=:positive, num_steps=n_16)
            A16             =       element2faces(E16)
            B16           .+=       Point{3, Float64}(0.0, -(0.5*TG+TW+0.5*TG), -0.5*Th_FL) 

            # dup_nodes       =       count_duplicate_nodes(B16)
            # println("Total duplicate nodes in B16 = $dup_nodes")

            ### CREATE SET - FL_LOAD / SUPPORT
            B16_SET_1       =       [p for p in B16 if p[2] == -(L + 0.5 * Gap_BtoB)]
            B16_SET_1_idx   =       findall(p -> p[2] == -(L + 0.5 * Gap_BtoB), B16)

            # =============================================================================
            # SCRIPT - END
            # =============================================================================
return B16, A16, E16, B16_SET_1, B16_SET_1_idx  
end

function comp_IBeam_v1(;
            TotalLength         = 197.0,                # Total length               - Beam to Beam
            Gap_BtoB            = 1.0,                  # Gap                        - Flange to Flange / Beam to Beam
            Edge_Support        = 50.0,                 # Edge Dist. to support      - Beam to Beam
            Th_FL               = 6.0,                  # Thickness                  - Flange
            Th_Web              = 6.0,                  # Thickness                  - Web
            nt_FL               = 2,                    # no. of teeth (each side)   - Flange 
            b_FL                = 50.0,                 # Width                      - Flange   
            Ht_Sec              = 50.0,                 # Height                     - Section
            TW                  = 10.0,                 # Tooth Width                - Flange
            TG                  = 11.0,                 # Tooth Gap                  - Flange  
            TD                  = 11.5,                 # Tooth Depth                - Flange
            r                   = 1.0,                  # Corner Radii               - Flange   
            g                   = 1.0,                  # Tolerance                  - b/w Tooth and Holes
            searchTol           = 1e-6,)                # Search Tolerance           - General
        # =============================================================================
        # SCRIPT - START
        # =============================================================================
        B16, A16, E16, B16_SET_1, B16_SET_1_idx = comp_FL_v1(TotalLength = TotalLength, Gap_BtoB = Gap_BtoB, Edge_Support = Edge_Support, Th_FL = Th_FL, Th_Web = Th_Web, nt_FL = nt_FL, b_FL = b_FL, 
        Ht_Sec = Ht_Sec, TW = TW, TG = TG, TD = TD, r = r, g = g, searchTol = searchTol,) 

        ### MESH CONTROL
        n_0             =       4                       
        n_1             =       4
        n_2             =       7                                                 # to link it
        n_3             =       5
        n_4             =       5
        n_5             =       n_0
        n_10            =       3
        n_11            =       n_1
        n_14            =       30
        n_16            =       4                                                # along Th_FL
        n_18            =       50  
        n_19            =       5   #   n_1+2                                            # along Th_Web

        ### DERIVED INPUTS
        Ht_Web          =       Ht_Sec-2*Th_FL         
        L               =       0.5*TotalLength-0.5*Gap_BtoB                    # Length of the Part excl. Gap_BtoB
        FL_B_L          =       L-(nt_FL*(TW + TG) + 0.5*TG - 0.5*Gap_BtoB)     # Bottom Length   - after the teeth section
        be_FL           =       b_FL - 2 * TD                                   # Effective Width - Flange
        nSections       =       nt_FL
        Pt_Load         =       TotalLength/6
        Pt_Support      =       L-Edge_Support
        searchTol       =       1e-6

        ### WEB    - LOFT # 19 20 
        B18_S1            =       filter(p -> isapprox(p[3],-0.5*Th_FL,atol=searchTol), B16)
        B18_S2            =       filter(p -> isapprox(p[1],-0.5*Th_Web,atol=searchTol), B18_S1)
        indSort           =       reverse(sortperm([v[2] for v in B18_S2]))
        B18_S2            =       B18_S2[indSort]
        B18_E2            =       [Point{3, Float64}(v[1], v[2], v[3]-Ht_Web) for v in B18_S2]
        A18, B18          =       loftlinear(B18_E2, B18_S2; num_steps = n_18, close_loop = false, face_type = :quad)

        E19, B19          =       extrudefaces(A18, B18; extent=Th_Web, direction=:negative, num_steps=n_19)
        A19               =       element2faces(E19)

        ### FL_BOT - LOFT # 21
        B21               =       [Point{3, Float64}(p[1], p[2], p[3] - (Ht_Web+Th_FL)) for p in B16]
        E21               =       copy(E16)
        A21               =       copy(A16)

        ### IBEAM  - LOFT # 22 (FL_TOP, WEB, FL_BOT)
        B22               =       [B16;B19;B21]
        B22             .+=       Point{3, Float64}(0.0, 0.0, +(0.5*Th_FL+0.5*Ht_Web))
        E22               =       vcat(E16, 
                                ISC.shift_elements(E19,length(B16)),
                                ISC.shift_elements(E21,length(B16)+length(B19)))
        A22               =       element2faces(E22)

        ### IBEAM - LOFT # 23 - REMOVING DUPLICATE NODES
        B23, E23, Old_New =      merge_nodes_and_update_elements(B22, E22)
        A23               =      element2faces(E23)

        # dup_nodes         =      count_duplicate_nodes(B23)
        # println("Total duplicate nodes in B23 = $dup_nodes")

        ### CREATE SET - Pt_Support
        B23_SET_1         =     [p for p in B23
                                if abs(p[2] - (-L + 50)) < 0.5*Edge_Support &&
                                abs(p[3] - (-0.5 * Ht_Sec)) <= searchTol]

        B23_SET_1_idx     =     findall(p -> abs(p[2] - (-L + 50)) < 0.5 * Edge_Support &&
                                abs(p[3] - (-0.5 * Ht_Sec)) <= searchTol, B23)

        ### CREATE SET - Pt_Load
        B23_SET_2         =     [p for p in B23
                                if abs(p[2] - (-Pt_Load)) < 0.5*Edge_Support &&
                                abs(p[3] - (0.5 * Ht_Sec)) <= searchTol]

        B23_SET_2_idx     =     findall(p -> abs(p[2] - (-Pt_Load)) < 0.5*Edge_Support &&
                                abs(p[3] - (0.5 * Ht_Sec)) <= searchTol, B23)

        ### CREATE SET - FL for Web Stiffener - R. Side
        B23_SET_3         =     [p for p in B23
                                if p[1] >= +(0.5*Th_Web) &&
                                (abs(p[3] - (-1*(0.5*Ht_Sec-Th_FL))) <= searchTol ||
                                abs(p[3]) <= (0.5*Ht_Sec - Th_FL) ||
                                abs(p[3] - (+1*(0.5*Ht_Sec-Th_FL))) <= searchTol)]

        B23_SET_3_idx     =     findall(p -> p[1] >= +(0.5*Th_Web) &&
                                (abs(p[3] - (-1*(0.5*Ht_Sec-Th_FL))) <= searchTol ||
                                abs(p[3]) <= (0.5*Ht_Sec - Th_FL) ||
                                abs(p[3] - (+1*(0.5*Ht_Sec-Th_FL))) <= searchTol),B23)

        ### CREATE SET - FL for Web Stiffener - L. Side
        B23_SET_4         =     [p for p in B23
                                if p[1] <= -(0.5*Th_Web) &&
                                (abs(p[3] - (-1*(0.5*Ht_Sec-Th_FL))) <= searchTol ||
                                abs(p[3]) <= (0.5*Ht_Sec - Th_FL) ||
                                abs(p[3] - (+1*(0.5*Ht_Sec-Th_FL))) <= searchTol)]

        B23_SET_4_idx   =       findall(p -> p[1] <= -(0.5*Th_Web) &&
                                (abs(p[3] - (-1*(0.5*Ht_Sec-Th_FL))) <= searchTol ||
                                abs(p[3]) <= (0.5*Ht_Sec - Th_FL) ||
                                abs(p[3] - (+1*(0.5*Ht_Sec-Th_FL))) <= searchTol),B23)

        # =============================================================================
        # SCRIPT - START
        # =============================================================================  
    return B22, A22
end

function comp_FL_Pantelis_v1(;
        TotalLength         =       500.0,    
        Gap_BtoB            =       0.0 ,                
        Edge_Support        =       50.0,                    
        Th_FL               =       9.6,        
        Th_Web              =       6.40,                                  
        b_FL                =       133.90,           
        Ht_Sec              =       206.80,                                
        TW                  =       30.0,        
        TG                  =       11.0,                              
        TD                  =       15.0,                    
        r                   =       1.0,                         
        g                   =       1.0,                       
        searchTol           =       1e-6,)
        # =============================================================================
        # SCRIPT - START
        # =============================================================================
        ### DERIVED INPUTS
        Ht_Web          =       Ht_Sec-2*Th_FL         
        L               =       0.5*TotalLength-0.5*Gap_BtoB                    # Length of the Part excl. Gap_BtoB
        FL_B_L          =       L-(nt_FL*(TW + TG) + 0.5*TG - 0.5*Gap_BtoB)     # Bottom Length   - after the teeth section
        be_FL           =       b_FL - 2 * TD                                   # Effective Width - Flange
        nSections       =       nt_FL
        Pt_Load         =       TotalLength/6
        Pt_Support      =       L-Edge_Support
        searchTol       =       1e-6

        ### MESH CONTROL
        n_0             =       4                       
        n_1             =       2
        n_2             =       5                                                 # to link it
        n_3             =       5
        n_4             =       2
        n_5             =       n_0
        n_10            =       3
        n_11            =       n_1
        n_14            =       30
        n_16            =       4                                                # along Th_FL
        n_18            =       50  
        n_19            =       1   #   n_1+2                                            # along Th_Web

        ### FL_TOP - LOFT # 1, EL # 1 - HI to ED           
        e0_h            =       0.5*TW+r
        e0_w            =       0.5*Th_Web
        V0_S            =       Vector{Point{3, Float64}}(undef, n_0)
        v_1             =       Point{3, Float64}(0.0, 0.0 , 0.0)
        v_2             =       Point{3, Float64}(0.0, e0_h , 0.0)
        for (i, w) in enumerate(range(0.0, 1.0, n_0))             
            V0_S[i]     =       (1.0 - w) * v_1 + w * v_2          
        end
        V0_E            =       [Point{3, Float64}(v[1]-e0_w, v[2], v[3]) for v in V0_S]
        A1, B1          =       loftlinear(V0_S, V0_E; num_steps = n_1, close_loop = false, face_type = :quad)
        ### FL_TOP - LOFT # 2, EL # 2 - HI to AB
        e1_h            =       0.5*TW+r
        e1_w            =       0.5*be_FL-0.5*TW-e0_w  
        V1_S            =       [Point{3, Float64}(v[1], v[2], v[3]) for v in V0_E]                               
        V1_E            =       [Point{3, Float64}(v[1]-e1_w, v[2], v[3]) for v in V1_S]
        A2, B2          =       loftlinear(V1_S, V1_E; num_steps = n_2, close_loop = false, face_type = :quad)

        ### FL_TOP - LOFT # 3, EL # 3 - AB BC to XY       
        e2_w            =       0.5*TW 
        Vr              =       [Point{3, Float64}(r * cosd(t), r * sind(t), 0.0) for t in range(0.0, 90, n_0+(n_0-1))]
        Vr              =       [Point{3, Float64}(v[1]-e0_w-e1_w-e1_h, v[2], v[3]) for v in Vr]
        Vb              =       Vector{Point{3, Float64}}(undef, n_0)  
        Vb1_E           =       [Point{3, Float64}(v[1], v[2], v[3]) for v in V1_E]
        Vb2_E           =       Vector{Point{3, Float64}}(undef, n_0)          
        v_1             =       Point{3, Float64}(-e0_w-e1_w        , e1_h , 0.0)
        v_2             =       Point{3, Float64}(-e0_w-e1_w-e1_h   , e1_h , 0.0)
        for (i, w) in enumerate(range(0.0, 1.0, n_0))              
            Vb2_E[i]    = (1.0 - w) * v_1 + w * v_2                
        end
        Vb              =       [Vb1_E;Vb2_E]
        Vb              =       unique(Vb, dims=1)
        A3, B3          =       loftlinear(Vb, Vr; num_steps = n_3, close_loop = false, face_type = :quad)

        ### FL_TOP - LOFT # 4, EL # 4 - HX to GF
        e3_h            =       0.5*TG-r   
        V3_S0           =       filter(p -> isapprox(p[2], 0.0 ,atol=searchTol), B1)
        V3_S1           =       filter(p -> isapprox(p[2], 0.0 ,atol=searchTol), B2)
        V3_S2           =       filter(p -> isapprox(p[2], 0.0 ,atol=searchTol), B3)
        V3_S            =       [V3_S0;V3_S1;V3_S2]
        V3_S            =       unique(V3_S, dims=1)
        V3_E            =       [Point{3, Float64}(v[1], v[2]-e3_h, v[3]) for v in V3_S]
        A4, B4          =       loftlinear(V3_S, V3_E; num_steps = n_4, close_loop = false, face_type = :quad)

        ### FL_TOP - LOFT # 5, EL # 5 - YC to DE
        e4_w            =       TD-r
        V4_S            =       filter(p -> isapprox(p[1], -e0_w-e1_w-e1_h ,atol=searchTol), B3)
        V4_E            =       [Point{3, Float64}(v[1]-e4_w, v[2], v[3]) for v in V4_S]
        A5, B5          =       loftlinear(V4_E, V4_S; num_steps = n_5, close_loop = false, face_type = :quad)

        ### FL_TOP - LOFT # 6 (1,2,3,4,5)
        B6              =       [B1; B2; B3; B4; B5]
        A6              =       [A1;
                                [f .+ length(B1) for f in A2];
                                [f .+ length(B1) .+ length(B2) for f in A3];
                                [f .+ length(B1) .+ length(B2) .+ length(B3) for f in A4];
                                [f .+ length(B1) .+ length(B2) .+ length(B3).+ length(B4) for f in A5]]
        A6, B6          =       mergevertices(A6, B6)
        # Eb5           =       boundaryedges(A6)                  
        B6            .+=       Point{3, Float64}(0.0, +e3_h, 0.0) # if i turn it off, why code does not work

        ### FL_TOP - LOFT # 7 (BASIC ELEMENT)
        B7              =       [B6;
                                [Point{3, Float64}(-v[1],  v[2], v[3]) for v in B6];
                                [Point{3, Float64}( v[1], -v[2], v[3]) for v in B6];
                                [Point{3, Float64}(-v[1], -v[2], v[3]) for v in B6]]    
        A7              =       [A6;
                                [reverse(f)  .+    length(B6)      for f in A6];
                                [reverse(f)  .+    length(B6)*2    for f in A6];
                                [f           .+    length(B6)*3    for f in A6]]

        ### FL_TOP - LOFT # 8 (n COPY OF BASIC ELEMENTS)                        
        nf              =       length(A7)
        nv              =       length(B7)
        A8              =       Vector{QuadFace{Int64}}(undef,nf*nSections)
        B8              =       Vector{Point{3,Float64}}(undef,nv*nSections)
        i_f             =       1
        i_v             =       1
        s               =       0
        for q in 1:nSections
            A8[i_f:i_f+nf-1]  = [f.+s for f in A7]
            B8[i_v: i_v+nv-1] = [Point{3, Float64}(v[1], v[2] - (q-1)*(TW + TG), v[3]) for v in B7]
            i_f += nf
            i_v += nv
            s   += nv
        end
        A8, B8          =       mergevertices(A8, B8)

        ### FL_TOP -  LOFT # 9 (TOP LEFT 1) - PMcG SCALED 
        B9              =       [[Point{3, Float64}(v[1],  -v[2]  .+ (TW + TG)-0.5*TG+r, v[3]) for v in B2];
                                [Point{3, Float64}(v[1],   -v[2]  .+ (TW + TG)-0.5*TG+r, v[3]) for v in B3];
                                [Point{3, Float64}(v[1],   -v[2]  .+ (TW + TG)-0.5*TG+r, v[3]) for v in B5]]

        A9              =       [[reverse(f) for f in A2];
                                [reverse(f) .+ length(B2) for f in A3]; 
                                [reverse(f) .+ length(B2) .+ length(B3) for f in A5]]
        A9, B9          =       mergevertices(A9, B9)

        ### FL_TOP -  LOFT # 10 (TOP LEFT 2) - PMcG SCALED 
        B10_S           =       filter(p -> isapprox(p[2], 0.5*TG + TW + r ,atol=searchTol), B9)
        indSort         =       reverse(sortperm([v[1] for v in B10_S]))
        B10_S           =       B10_S[indSort]
        B10_E           =       [Point{3, Float64}(v[1],  v[2]  .+ (0.5*TG - 0.5*Gap_BtoB - r), v[3]) for v in B10_S]
        indSort         =       reverse(sortperm([v[1] for v in B10_S]))
        B10_E           =       B10_E[indSort]
        A10,B10         =       loftlinear(B10_E, B10_S ; num_steps = n_10, close_loop = false, face_type = :quad)

        ### FL_TOP -  LOFT # 11 (TOP LEFT 3) - PMcG SCALED
        B11_S1          =       filter(p -> isapprox(p[1], -e0_w ,atol=searchTol), B9)
        B11_S2          =       filter(p -> isapprox(p[1], -e0_w ,atol=searchTol), B10)
        B11_S           =       [B11_S1;B11_S2]
        indSort         =       reverse(sortperm([v[2] for v in B11_S]))
        B11_S           =       B11_S[indSort]
        B11_S           =       unique(B11_S, dims=1)
        B11_E           =       [Point{3, Float64}(v[1]+e0_w, v[2], v[3]) for v in B11_S]
        indSort         =       reverse(sortperm([v[2] for v in B11_E]))
        B11_E           =       B11_E[indSort]
        A11,B11         =       loftlinear(B11_S, B11_E; num_steps = n_11, close_loop = false, face_type = :quad)

        ### FL_TOP -  LOFT # 12 (9 10 11) - PMcG SCALED
        B12             =       [B9; B10; B11]                    
        A12             =       [A9;
                                [f .+ length(B9) for f in A10];
                                [f .+ length(B9) .+ length(B10) for f in A11]]
        A12, B12        =       mergevertices(A12, B12)
        # Eb11          =       boundaryedges(A12)  

        ### FL_TOP -  LOFT # 13 (12 12) - PMcG SCALED
        B13             =       [B12;
                                [Point{3, Float64}(-v[1],  v[2], v[3]) for v in B12]]     # top right    # bottom left
        A13             =       [A12;
                                [reverse(f)  .+    length(B12)     for f in A12]]
        A13, B13        =       mergevertices(A13, B13)
        # Eb_N2         =       boundaryedges(A13)

        ### FL_TOP - LOFT # 14 - BOTTOM PART
        B14_S           =       filter(p -> isapprox(p[2],-(nSections-0.5)*(TG+TW),atol=searchTol), B8)
        indSort         =       reverse(sortperm([v[1] for v in B14_S]))
        B14_S           =       B14_S[indSort]
        B14_E           =       [Point{3, Float64}(v[1], -(FL_B_L-0.5*TW)+v[2] , v[3]) for v in B14_S]
        A14, B14        =       loftlinear(B14_S, B14_E; num_steps = n_14, close_loop = false, face_type = :quad)
        ### FL_TOP - LOFT # 15_1 
        B15_1               =       [B8; B14]
        A15_1               =       [A8;
                                    [f .+ length(B8) for f in A14]]
        A15_1, B15_1        =       mergevertices(A15_1, B15_1) 

        ### FL_TOP - LOFT # 15_2
        B15_2_S             =       filter(p -> isapprox(
                                    p[2], 0.5*TG+0.5*TW, atol = searchTol)      && 
                                    p[1] > ((0.5*b_FL) - (b_FL/3.82) -(0.5*TW)) && 
                                    p[1] < ((0.5*b_FL) - (b_FL/3.82)) ,B15_1)

        B15_2_E1            =       [Point{3, Float64}(p[1], p[2] - TD, p[3]) for p in B15_2_S]
        B15_2_E1_idx        =       unique([nearest_node_id(p, B15_1) for p in B15_2_E1])
        B15_2_E             =       B15_1[B15_2_E1_idx]
                    
        B15_2               =       filter(p ->
                                    p[1] >= minimum(p[1] for p in B15_2_S) - searchTol &&
                                    p[1] <= maximum(p[1] for p in B15_2_S) + searchTol &&
                                    p[2] >= minimum(p[2] for p in vcat(B15_2_S, B15_2_E)) - searchTol &&
                                    p[2] <= maximum(p[2] for p in vcat(B15_2_S, B15_2_E)) + searchTol,B15_1)

        B15_3               =       [Point{3, Float64}(p[1], -p[2] + (TG + TW), p[3]) for p in B15_2]
        A15_4, B15_4        =       quad_faces_nodes(B15_3; dir1 = 1, dir2 = 2, tol = searchTol)

        ### FL_TOP - LOFT # 15_5 (Added Pantelis 1 Teeth)
        B15_5               =       [B15_1; B15_4]
        A15_5               =       [A15_1;
                                    [f .+ length(B15_1) for f in A15_4]]
        A15_5, B15_5        =       mergevertices(A15_5, B15_5)

        ### FL_TOP - LOFT # 15_6 
        B15_6               =       [Point{3, Float64}(-p[1], p[2], p[3]) for p in B15_2]  

        ### FL_TOP - LOFT # 16
        A15_Pan, B15_Pan    =       remove_faces_inside_nodes(A15_5, B15_5, B15_6; searchTol = searchTol)

        ### FL_TOP - LOFT # 16 - FINAL EXTRUDE
        E16_Pan,B16_Pan     =       extrudefaces(A15_Pan, B15_Pan; extent=Th_FL, direction=:positive, num_steps=n_16)
        A16_Pan             =       element2faces(E16_Pan)
        B16_Pan           .+=       Point{3, Float64}(0.0, -(0.5*TG+0.5*TW+0.5*Gap_BtoB), -0.5*Th_FL) 

        # =============================================================================
        # SCRIPT - END
        # =============================================================================
return B16_Pan, A16_Pan, E16_Pan
end