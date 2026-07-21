# FUNCTIONS - CREATE GEOMETRY
function quad_faces_nodes(V; dir1 = 1, dir2 = 2, tol = 1e-6)
    # Sort nodes by second direction, then first direction
    V_sorted = sort(V, by = p -> (round(p[dir2] / tol), round(p[dir1] / tol)))

    xvals    = unique(round.(getindex.(V_sorted, dir1) ./ tol))
    yvals    = unique(round.(getindex.(V_sorted, dir2) ./ tol))

    nx       = length(xvals)
    ny       = length(yvals)

    if length(V_sorted) != nx * ny
        error("V does not form a complete rectangular grid")
    end

    F = QuadFace{Int64}[]

    for j in 1:ny-1
        for i in 1:nx-1
            n1 = (j - 1) * nx + i
            n2 = n1 + 1
            n3 = n1 + nx + 1
            n4 = n1 + nx

            push!(F, QuadFace{Int64}(n1, n2, n3, n4))
        end
    end

    return F, V_sorted
end

function nearest_node_id(p, V)
    distances = [norm(v - p) for v in V]
    return argmin(distances)
end

function remove_faces_inside_nodes(A_original, B_original, B_remove; searchTol = 1e-6)
    function face_center(f, V)
        ids = Int.(collect(f))
        pts = V[ids]

        Point{3, Float64}(
            mean(p[1] for p in pts),
            mean(p[2] for p in pts),
            mean(p[3] for p in pts)
        )
    end

    x_min = minimum(p[1] for p in B_remove)
    x_max = maximum(p[1] for p in B_remove)

    y_min = minimum(p[2] for p in B_remove)
    y_max = maximum(p[2] for p in B_remove)

    z_ref = mean(p[3] for p in B_remove)

    A_keep_old = filter(f -> begin
            c = face_center(f, B_original)

            inside =
                c[1] >= x_min - searchTol &&
                c[1] <= x_max + searchTol &&
                c[2] >= y_min - searchTol &&
                c[2] <= y_max + searchTol &&
                isapprox(c[3], z_ref, atol = searchTol)

            !inside
        end,
        A_original
    )

    if isempty(A_keep_old)
        error("All faces were removed. Check selected removal nodes or searchTol.")
    end

    B_keep_ids = sort(unique(vcat([Int.(collect(f)) for f in A_keep_old]...)))

    B_new = B_original[B_keep_ids]

    old_to_new = Dict(old_id => new_id for (new_id, old_id) in enumerate(B_keep_ids))

    A_new = [
        QuadFace{Int64}(
            old_to_new[Int(collect(f)[1])],
            old_to_new[Int(collect(f)[2])],
            old_to_new[Int(collect(f)[3])],
            old_to_new[Int(collect(f)[4])]
        )
        for f in A_keep_old
    ]

    return A_new, B_new
end

function merge_nodes_and_update_elements(V, E; tol=1e-6)
    node_dict = Dict{Tuple{Int,Int,Int}, Int}()
    Vnew = Point{3,Float64}[]
    old_to_new = Vector{Int}(undef, length(V))

    for (i, p) in enumerate(V)
        key = (
            round(Int, p[1]/tol),
            round(Int, p[2]/tol),
            round(Int, p[3]/tol)
        )

        if haskey(node_dict, key)
            old_to_new[i] = node_dict[key]
        else
            new_id = length(Vnew) + 1
            push!(Vnew, p)
            node_dict[key] = new_id
            old_to_new[i] = new_id
        end
    end

    Enew = [
        Hex8{Int64}(
            old_to_new[e[1]], old_to_new[e[2]], old_to_new[e[3]], old_to_new[e[4]],
            old_to_new[e[5]], old_to_new[e[6]], old_to_new[e[7]], old_to_new[e[8]]
        )
        for e in E
    ]

    return Vnew, Enew, old_to_new
end

function count_duplicate_nodes(V; tol=1e-6)
    coord_count = Dict{Tuple{Int,Int,Int}, Int}()
    for p in V
        key = (
            round(Int, p[1]/tol),
            round(Int, p[2]/tol),
            round(Int, p[3]/tol)
        )
        coord_count[key] = get(coord_count, key, 0) + 1
    end
    duplicate_nodes = sum(v -> v - 1, values(coord_count))
    return duplicate_nodes
end

#cut_hole : Snap Rectangle 1 to existing mesh lines. Remove internal nodes and faces. Create inward-offset Rectangle 2. Create the circular boundary. Loft the two transition regions. Merge everything once at the end.
function cut_hole_yz(
    A1,
    B1;
    center,
    radius,
    effect_scale=6.0,
    offset_scale=1.5,
    n_radial=5,
    tol=1e-6)

    xc1, yc1, zc1 = center
    r = radius

    ### RECTANGLE 1: IMAGINARY

    Effect      = effect_scale * r
    Rec1_w      = Effect
    Rec1_h      = Effect
    offset_Rec2 = offset_scale * r

    w1_Rec1 = yc1 - 0.5 * Rec1_w
    w2_Rec1 = yc1 + 0.5 * Rec1_w
    h1_Rec1 = zc1 - 0.5 * Rec1_h
    h2_Rec1 = zc1 + 0.5 * Rec1_h

    Pt_TL = Point{3, Float64}(xc1, w1_Rec1, h2_Rec1)
    Pt_BL = Point{3, Float64}(xc1, w1_Rec1, h1_Rec1)
    Pt_TR = Point{3, Float64}(xc1, w2_Rec1, h2_Rec1)
    Pt_BR = Point{3, Float64}(xc1, w2_Rec1, h1_Rec1)

    ### SNAP RECTANGLE 1 CORNERS TO EXISTING MESH NODES

    closest_B1_node(p) = B1[
        argmin([
            (q[1] - p[1])^2 +
            (q[2] - p[2])^2 +
            (q[3] - p[3])^2
            for q in B1
        ])
    ]

    Pt_TL1 = closest_B1_node(Pt_TL)
    Pt_BL1 = closest_B1_node(Pt_BL)
    Pt_TR1 = closest_B1_node(Pt_TR)
    Pt_BR1 = closest_B1_node(Pt_BR)

    w1_Rec1 = Pt_TL1[2]
    w2_Rec1 = Pt_TR1[2]
    h1_Rec1 = Pt_BL1[3]
    h2_Rec1 = Pt_TL1[3]

    ### RECTANGLE 1 EDGES: EXISTING MESH NODES

    Edge_Rec1_L = sort(
        [p for p in B1 if
         isapprox(p[2], w1_Rec1; atol=tol) &&
         h1_Rec1-tol <= p[3] <= h2_Rec1+tol],
        by=p -> -p[3]
    )

    Edge_Rec1_R = sort(
        [p for p in B1 if
         isapprox(p[2], w2_Rec1; atol=tol) &&
         h1_Rec1-tol <= p[3] <= h2_Rec1+tol],
        by=p -> p[3]
    )

    Edge_Rec1_B = sort(
        [p for p in B1 if
         isapprox(p[3], h1_Rec1; atol=tol) &&
         w1_Rec1-tol <= p[2] <= w2_Rec1+tol],
        by=p -> p[2]
    )

    Edge_Rec1_T = sort(
        [p for p in B1 if
         isapprox(p[3], h2_Rec1; atol=tol) &&
         w1_Rec1-tol <= p[2] <= w2_Rec1+tol],
        by=p -> -p[2]
    )

    @assert minimum(length.(
        (Edge_Rec1_L, Edge_Rec1_R, Edge_Rec1_B, Edge_Rec1_T)
    )) >= 2 "Not enough existing mesh nodes around the selected rectangle."

    Edge_Rec1 = [
        Edge_Rec1_B[1:end-1];
        Edge_Rec1_R[1:end-1];
        Edge_Rec1_T[1:end-1];
        Edge_Rec1_L
    ]

    ### REMOVE INTERNAL NODES AND FACES

    inside_rec1(p) =
        w1_Rec1+tol < p[2] < w2_Rec1-tol &&
        h1_Rec1+tol < p[3] < h2_Rec1-tol

    keep_node = [!inside_rec1(p) for p in B1]

    old_to_new = zeros(Int, length(B1))
    old_to_new[keep_node] = 1:count(keep_node)

    B1_blank = B1[keep_node]
    A1_blank = eltype(A1)[]

    for f in A1
        ids = collect(Tuple(f))

        if all(keep_node[ids])
            new_ids = old_to_new[ids]
            push!(A1_blank, typeof(f)(new_ids...))
        end
    end

    A1_blank, B1_blank = mergevertices(A1_blank, B1_blank)

    ### RECTANGLE 2: INWARD OFFSET

    w1_Rec2 = w1_Rec1 + offset_Rec2
    w2_Rec2 = w2_Rec1 - offset_Rec2
    h1_Rec2 = h1_Rec1 + offset_Rec2
    h2_Rec2 = h2_Rec1 - offset_Rec2

    @assert w1_Rec2 < w2_Rec2 "Offset is too large in the Y direction."
    @assert h1_Rec2 < h2_Rec2 "Offset is too large in the Z direction."

    Pt_TL2 = Point{3, Float64}(xc1, w1_Rec2, h2_Rec2)
    Pt_BL2 = Point{3, Float64}(xc1, w1_Rec2, h1_Rec2)
    Pt_TR2 = Point{3, Float64}(xc1, w2_Rec2, h2_Rec2)
    Pt_BR2 = Point{3, Float64}(xc1, w2_Rec2, h1_Rec2)

    Edge_Rec2_L = [
        (1.0-w)*Pt_TL2 + w*Pt_BL2
        for w in range(0.0, 1.0, length=length(Edge_Rec1_L))
    ]

    Edge_Rec2_R = [
        (1.0-w)*Pt_BR2 + w*Pt_TR2
        for w in range(0.0, 1.0, length=length(Edge_Rec1_R))
    ]

    Edge_Rec2_B = [
        (1.0-w)*Pt_BL2 + w*Pt_BR2
        for w in range(0.0, 1.0, length=length(Edge_Rec1_B))
    ]

    Edge_Rec2_T = [
        (1.0-w)*Pt_TR2 + w*Pt_TL2
        for w in range(0.0, 1.0, length=length(Edge_Rec1_T))
    ]

    Edge_Rec2 = [
        Edge_Rec2_B[1:end-1];
        Edge_Rec2_R[1:end-1];
        Edge_Rec2_T[1:end-1];
        Edge_Rec2_L
    ]

    ### CIRCULAR ARCS

    Arc_Cir_L = [
        Point{3, Float64}(xc1, yc1+r*cos(t), zc1+r*sin(t))
        for t in range(0.75pi, 1.25pi, length=length(Edge_Rec2_L))
    ]

    Arc_Cir_R = [
        Point{3, Float64}(xc1, yc1+r*cos(t), zc1+r*sin(t))
        for t in range(1.75pi, 2.25pi, length=length(Edge_Rec2_R))
    ]

    Arc_Cir_B = [
        Point{3, Float64}(xc1, yc1+r*cos(t), zc1+r*sin(t))
        for t in range(1.25pi, 1.75pi, length=length(Edge_Rec2_B))
    ]

    Arc_Cir_T = [
        Point{3, Float64}(xc1, yc1+r*cos(t), zc1+r*sin(t))
        for t in range(0.25pi, 0.75pi, length=length(Edge_Rec2_T))
    ]

    Arc_Cir = [
        Arc_Cir_B[1:end-1];
        Arc_Cir_R[1:end-1];
        Arc_Cir_T[1:end-1];
        Arc_Cir_L
    ]

    ### LOFT TRANSITIONS

    A_Rec1, B_Rec1 = loftlinear(
        Edge_Rec1,
        Edge_Rec2;
        num_steps=n_radial,
        close_loop=false,
        face_type=:quad
    )

    A_R_Arc, B_R_Arc = loftlinear(
        Edge_Rec2,
        Arc_Cir;
        num_steps=n_radial,
        close_loop=false,
        face_type=:quad
    )

    ### COMBINE AND MERGE

    shift_faces(F, shift) = [
        typeof(f)((collect(Tuple(f)) .+ shift)...)
        for f in F
    ]

    B_combined = [B1_blank; B_Rec1; B_R_Arc]

    A_combined = [
        A1_blank;
        shift_faces(A_Rec1, length(B1_blank));
        shift_faces(A_R_Arc, length(B1_blank) + length(B_Rec1))
    ]

    A_final, B_final = mergevertices(A_combined, B_combined)

    boundaries = (
        rec1=Edge_Rec1,
        rec2=Edge_Rec2,
        circle=Arc_Cir
    )

    return A_final, B_final, boundaries
end # TO DO - for XY and YZ

function create_rect(;
    wd = 100.0,
    ht = 200.0,
    n_wd = 20,
    n_ht = 40,
    center    = (0.0, 0.0, 0.0),
    plane = :XY,)

    # INPUTS - DERIVED
    center = Point{3, Float64}(center...)
    cx, cy, cz = center

    # CREATE - POINTS
    if     plane == :YZ 
        p1 = Point{3, Float64}(cx, cy - wd / 2, cz - ht / 2)
        p2 = Point{3, Float64}(cx, cy + wd / 2, cz - ht / 2)
    elseif plane == :XY
        p1 = Point{3, Float64}(cx - wd / 2, cy - ht / 2, cz)
        p2 = Point{3, Float64}(cx + wd / 2, cy - ht / 2, cz)
    else
        error("Unsupported plane: $plane")
    end

    # CREATE - EDGES
    edge_start = [(1.0 - w) * p1 + w * p2
                for w in range(0.0, 1.0, length = n_wd + 1)]

    if plane == :YZ
    edge_end = [Point{3,Float64}(p[1], p[2], p[3] + ht) for p in edge_start]
    elseif plane == :XY
    edge_end = [Point{3,Float64}(p[1], p[2] + ht, p[3]) for p in edge_start]
    end

    # CREATE - LOFT
    faces, vertices = loftlinear(edge_end, edge_start; num_steps = n_ht, close_loop = false, face_type = :quad)

    return faces, vertices
end

function remove_duplicate_nodes(V; digits = 8)
    key(p) = (
        round(p[1], digits = digits),
        round(p[2], digits = digits),
        round(p[3], digits = digits),
    )

    seen = Dict{Tuple{Float64, Float64, Float64}, Int}()
    V_unique = Point{3, Float64}[]
    old_to_new = Vector{Int}(undef, length(V))

    for (i, p) in enumerate(V)
        k = key(p)

        if haskey(seen, k)
            old_to_new[i] = seen[k]
        else
            push!(V_unique, p)
            seen[k] = length(V_unique)
            old_to_new[i] = length(V_unique)
        end
    end

    return V_unique, old_to_new
end

function create_circle(; 
    Center = (0.0 , 0.0 , 0.0),
    r      = 50.0, 
    n_Edge1= 6,
    )

    Center = Point{3, Float64}(Center...)

    # INPUTS - DERIVED
    n_Edge2         =       ceil(Int, 0.5*n_Edge1)
    fac1            =       0.50         
    fac2            =       0.75

    # CREATE - IMAGINARY RECTANGLE


    Point1          =       Point{3, Float64}(0.0        , fac1*r     ,0.0)
    Point2          =       Point{3, Float64}(fac2*fac1*r, fac1*r    ,0.0)
    Point3          =       Point{3, Float64}(fac1*r     , fac2*fac1*r,0.0)
    Point4          =       Point{3, Float64}(fac1*r     , 0.0         ,0.0)

    Edge1           =       Vector{Point{3, Float64}}(undef, n_Edge1)
    for (i, w) in enumerate(range(0.0, 1.0, n_Edge1))
        Edge1[i]    =       (1.0 - w) * Point1 + w * Point2
    end

    Edge2           =        Vector{Point{3, Float64}}(undef, n_Edge2)
    for (i, w) in enumerate(range(0.0, 1.0, n_Edge2))
        Edge2[i]    =       (1.0 - w) * Point2 + w * Point3
    end

    Edge12          =       [Edge1 ; Edge2]
    Edge12          =       remove_duplicate_nodes(Edge12; digits = 8)[1]

    indSort         =       reverse(sortperm([v[1] for v in Edge12]))
    Edge12          =       Edge12[indSort]

    Edge4           =       Vector{Point{3, Float64}}(undef, length(Edge12))
    for (i, w) in enumerate(range(0.0, 1.0, length(Edge12)))
        Edge4[i]    =       (1.0 - w) * Point4 + w * Center
    end

    # CREATE - CIRCLE - QUARTER - P1
    F1, V1          =       loftlinear(Edge12, Edge4 ; num_steps = n_Edge1, close_loop = false, face_type = :quad)
    n_edge          =       length(Edge4)                   # to understand this and check if it could be added in Commodo.
    n_rows          =       div(length(V1), n_edge)
    Vgrid           =       reshape(V1, n_edge, n_rows)
    Edge3           =       collect(Vgrid[1, :])
    Edge5           =       collect(Vgrid[end, :])

    Edge123         =       [Edge12; Edge3]

    Edge123         =       remove_duplicate_nodes(Edge123; digits = 8)[1]

    indSort         =       reverse(sortperm([v[1] for v in Edge123]))
    Edge123         =       Edge123[indSort]

    indSort         =       reverse(sortperm([v[2] for v in Edge123]))
    Edge123         =       Edge123[indSort]

    # CREATE - CIRCLE - QUARTER - P2
    Br              =      [Point{3, Float64}(r * sind(t), r * cosd(t), 0.0)
                            for t in range(0.0, 90, length(Edge123))]

    F2, V2          =       loftlinear(Edge123, Br ; num_steps = n_Edge1, close_loop = false, face_type = :quad)
    F2, V2          =       mergevertices(F2, V2)

    # CREATE - CIRCLE - QUARTER - P1 + P2
    V               =       [V1;V2]
    F               =       [F1; 
                            [(f)  .+    length(V1)      for f in F2]]

    F, V            =       mergevertices(F, V)

    # CREATE - CIRCLE

    V_Cir1          =       [V;
                            [Point{3, Float64}(0.0-v[1],     v[2], v[3]) for v in V];
                            [Point{3, Float64}(0.0-v[1], 0.0-v[2], v[3]) for v in V];
                            [Point{3, Float64}(    v[1], 0.0-v[2], v[3]) for v in V]]    

    F_Cir1          =       [F;
                            [reverse(f)  .+    length(V)      for f in F];
                            [       (f)  .+    length(V)*2    for f in F];
                            [reverse(f)  .+    length(V)*3    for f in F]]

    F_Cir1, V_Cir1  =        mergevertices(F_Cir1, V_Cir1)

    return (F_Cir1, V_Cir1)
end # TO DO - for XY and YZ

function cut_rect(A, B;
    Center = (0.0, -30.0, 0.0),
    Rec1_w = 40.0,
    Rec1_h = 40.0,
    offset = 20.0,
    n_off = 4,
    plane = :XY,
    Edge = :L,
    tol = 1e-6,)

    # INPUTS - DERIVED
    xc1, yc1, zc1       = Center
 
    # REC 1 - EXACT REQUIRED
    if     plane == :YZ
        w1_Rec1             =   yc1-0.5*Rec1_w
        w2_Rec1             =   yc1+0.5*Rec1_w
        h1_Rec1             =   zc1-0.5*Rec1_h
        h2_Rec1             =   zc1+0.5*Rec1_h
    elseif plane == :XY
        w1_Rec1             =   xc1-0.5*Rec1_w
        w2_Rec1             =   xc1+0.5*Rec1_w
        h1_Rec1             =   yc1-0.5*Rec1_h
        h2_Rec1             =   yc1+0.5*Rec1_h
    else
    error("Unsupported plane: $plane")
    end

    if     plane == :YZ 
        Pt_TL1              =   Point{3 , Float64}(xc1, w1_Rec1, h2_Rec1)
        Pt_BL1              =   Point{3 , Float64}(xc1, w1_Rec1, h1_Rec1)
        Pt_TR1              =   Point{3 , Float64}(xc1, w2_Rec1, h2_Rec1)
        Pt_BR1              =   Point{3 , Float64}(xc1, w2_Rec1, h1_Rec1)
    elseif plane == :XY
        Pt_TL1              =   Point{3 , Float64}(w1_Rec1, h2_Rec1, zc1)
        Pt_BL1              =   Point{3 , Float64}(w1_Rec1, h1_Rec1, zc1)
        Pt_TR1              =   Point{3 , Float64}(w2_Rec1, h2_Rec1, zc1)
        Pt_BR1              =   Point{3 , Float64}(w2_Rec1, h1_Rec1, zc1)
    else
        error("Unsupported plane: $plane")
    end

    # REC 2 - IMAGINARY
    w1_Rec2             =   w1_Rec1-0.5*offset
    w2_Rec2             =   w2_Rec1+0.5*offset
    h1_Rec2             =   h1_Rec1-0.5*offset
    h2_Rec2             =   h2_Rec1+0.5*offset

    if     plane == :YZ 
        Pt_TL2              =   Point{3 , Float64}(xc1, w1_Rec2, h2_Rec2)
        Pt_BL2              =   Point{3 , Float64}(xc1, w1_Rec2, h1_Rec2)
        Pt_TR2              =   Point{3 , Float64}(xc1, w2_Rec2, h2_Rec2)
        Pt_BR2              =   Point{3 , Float64}(xc1, w2_Rec2, h1_Rec2)
    elseif plane == :XY
        Pt_TL2              =   Point{3 , Float64}(w1_Rec2, h2_Rec2, zc1)
        Pt_BL2              =   Point{3 , Float64}(w1_Rec2, h1_Rec2, zc1)
        Pt_TR2              =   Point{3 , Float64}(w2_Rec2, h2_Rec2, zc1)
        Pt_BR2              =   Point{3 , Float64}(w2_Rec2, h1_Rec2, zc1)
    else
        error("Unsupported plane: $plane")
    end

    # REC 3 - CORNERS - SNAP EXIST MESH CORNERS (ACTUALLY = REC 2)
    xs                      =   sort(unique(p[1] for p in B))
    ys                      =   sort(unique(p[2] for p in B))
    zs                      =   sort(unique(p[3] for p in B))
    if     plane == :YZ 
        w1_Rec3             =   ys[argmin(abs.(ys .- w1_Rec2))]
        w2_Rec3             =   ys[argmin(abs.(ys .- w2_Rec2))]
        h1_Rec3             =   zs[argmin(abs.(zs .- h1_Rec2))]
        h2_Rec3             =   zs[argmin(abs.(zs .- h2_Rec2))]
    elseif plane == :XY
        w1_Rec3             =   xs[argmin(abs.(xs .- w1_Rec2))]
        w2_Rec3             =   xs[argmin(abs.(xs .- w2_Rec2))]
        h1_Rec3             =   ys[argmin(abs.(ys .- h1_Rec2))]
        h2_Rec3             =   ys[argmin(abs.(ys .- h2_Rec2))]
    else
        error("Unsupported plane: $plane")
    end

    if     plane == :YZ 
    Pt_TL3              =   Point{3 , Float64}(xc1, w1_Rec3, h2_Rec3)
    Pt_BL3              =   Point{3 , Float64}(xc1, w1_Rec3, h1_Rec3)
    Pt_TR3              =   Point{3 , Float64}(xc1, w2_Rec3, h2_Rec3)
    Pt_BR3              =   Point{3 , Float64}(xc1, w2_Rec3, h1_Rec3)
    elseif plane == :XY
    Pt_TL3              =   Point{3 , Float64}(w1_Rec3, h2_Rec3, zc1)
    Pt_BL3              =   Point{3 , Float64}(w1_Rec3, h1_Rec3, zc1)
    Pt_TR3              =   Point{3 , Float64}(w2_Rec3, h2_Rec3, zc1)
    Pt_BR3              =   Point{3 , Float64}(w2_Rec3, h1_Rec3, zc1)
    else
        error("Unsupported plane: $plane")
    end

    if     plane == :YZ 
    Edge_Rec3_L         =   sort([p for p in B if isapprox(p[2], w1_Rec3; atol=tol) &&         #TOP    -    BOT
                            h1_Rec3 - tol <= p[3] <= h2_Rec3 + tol], by = p -> -p[3])
    Edge_Rec3_R         =   sort([p for p in B if isapprox(p[2], w2_Rec3; atol=tol) &&         #BOT    -    TOP
                            h1_Rec3 - tol <= p[3] <= h2_Rec3 + tol], by = p -> +p[3])
    Edge_Rec3_B         =   sort([p for p in B if isapprox(p[3], h1_Rec3; atol=tol) &&         #LEFT   -    RIGHT
                            w1_Rec3 - tol <= p[2] <= w2_Rec3 + tol],by = p -> p[2])
    Edge_Rec3_T         =   sort([p for p in B if isapprox(p[3], h2_Rec3; atol=tol) &&         #RIGHT  -    LEFT
                            w1_Rec3 - tol <= p[2] <= w2_Rec3 + tol], by = p -> -p[2])
    elseif plane == :XY
    Edge_Rec3_L         =   sort([p for p in B if isapprox(p[1], w1_Rec3; atol=tol) && 
                            h1_Rec3-tol <= p[2] <= h2_Rec3+tol], by=p -> -p[2])
    Edge_Rec3_R         =   sort([p for p in B if isapprox(p[1], w2_Rec3; atol=tol) && 
                            h1_Rec3-tol <= p[2] <= h2_Rec3+tol], by=p ->  p[2])
    Edge_Rec3_B         =   sort([p for p in B if isapprox(p[2], h1_Rec3; atol=tol) && 
                            w1_Rec3-tol <= p[1] <= w2_Rec3+tol], by=p ->  p[1])
    Edge_Rec3_T         =   sort([p for p in B if isapprox(p[2], h2_Rec3; atol=tol) && 
                            w1_Rec3-tol <= p[1] <= w2_Rec3+tol], by=p -> -p[1])
    else
        error("Unsupported plane: $plane")
    end

    # REMOVE EVERYTHING INSIDE THE REC 3

    inside_rec1(p)      = plane == :YZ ? (w1_Rec3+tol < p[2] < w2_Rec3-tol && h1_Rec3+tol< p[3] < h2_Rec3-tol) : 
                        plane == :XY ? (w1_Rec3+tol < p[1] < w2_Rec3-tol && h1_Rec3+tol < p[2] < h2_Rec3-tol) : 
                        error("Unsupported plane: $plane")

    keep_node           =   [!inside_rec1(p) for p in B]

    # Retain only faces whose nodes are all outside the rectangle
    kept_faces = [
        f for f in A
        if all(keep_node[Int(i)] for i in f)
    ]

    # Determine which original nodes are actually used by retained faces
    used_old = sort(unique(collect(Iterators.flatten(kept_faces))))

    # Map original node indices to compact indices
    old_to_new = zeros(Int, length(B))
    old_to_new[used_old] = collect(eachindex(used_old))

    # Keep only nodes referenced by retained faces
    B_blank = B[used_old]

    # Rebuild retained faces using compact node indices
    A_blank = [
        typeof(f)((old_to_new[Int(i)] for i in f)...)
        for f in kept_faces
    ]

    # Merge any coincident vertices
    A_blank, B_blank = mergevertices(A_blank, B_blank)

    # REC 1 - EXACT REQUIRED
    Edge_Rec1_L         =  [(1.0-w)*Pt_TL1 + w*Pt_BL1 for w in range(0.0, 1.0, length=length(Edge_Rec3_L))]  #TOP    -    BOT
    Edge_Rec1_R         =  [(1.0-w)*Pt_BR1 + w*Pt_TR1 for w in range(0.0, 1.0, length=length(Edge_Rec3_R))]  #BOT    -    TOP
    Edge_Rec1_B         =  [(1.0-w)*Pt_BL1 + w*Pt_BR1 for w in range(0.0, 1.0, length=length(Edge_Rec3_B))]  #LEFT   -    RIGHT
    Edge_Rec1_T         =  [(1.0-w)*Pt_TR1 + w*Pt_TL1 for w in range(0.0, 1.0, length=length(Edge_Rec3_T))]  #RIGHT  -    LEFT

    # CUT REC - EDGE - I 
    Edge_Rec3           =  [Edge_Rec3_B[1:end-1]; Edge_Rec3_R[1:end-1]; Edge_Rec3_T[1:end-1]; Edge_Rec3_L[1:end]]
    Edge_Rec1           =  [Edge_Rec1_B[1:end-1]; Edge_Rec1_R[1:end-1]; Edge_Rec1_T[1:end-1]; Edge_Rec1_L[1:end]]
    A_Rec1_I, B_Rec1_I  =  loftlinear(Edge_Rec1, Edge_Rec3 ; num_steps=n_off, close_loop=false, face_type=:quad)

    # CUT REC - EDGE - R 
    Edge_Rec3           =   [Edge_Rec3_T[1:end-1]; Edge_Rec3_L[1:end-1]; Edge_Rec3_B[1:end]]
    Edge_Rec1           =   [Edge_Rec1_T[1:end-1]; Edge_Rec1_L[1:end-1]; Edge_Rec1_B[1:end]]
    A_Rec1_R, B_Rec1_R  =   loftlinear(Edge_Rec1, Edge_Rec3 ; num_steps=n_off, close_loop=false, face_type=:quad)

    # CUT REC - EDGE - L
    Edge_Rec3           =   [Edge_Rec3_B[1:end-1]; Edge_Rec3_R[1:end-1]; Edge_Rec3_T[1:end]]
    Edge_Rec1           =   [Edge_Rec1_B[1:end-1]; Edge_Rec1_R[1:end-1]; Edge_Rec1_T[1:end]]
    A_Rec1_L, B_Rec1_L  =   loftlinear(Edge_Rec1, Edge_Rec3 ; num_steps=n_off, close_loop=false, face_type=:quad)

    # CUT REC - EDGE - T
    Edge_Rec3           =   [Edge_Rec3_L[1:end-1]; Edge_Rec3_B[1:end-1]; Edge_Rec3_R[1:end]]
    Edge_Rec1           =   [Edge_Rec1_L[1:end-1]; Edge_Rec1_B[1:end-1]; Edge_Rec1_R[1:end]]
    A_Rec1_T, B_Rec1_T  =   loftlinear(Edge_Rec1, Edge_Rec3 ; num_steps=n_off, close_loop=false, face_type=:quad)

    # CUT REC - EDGE - B
    Edge_Rec3           =   [Edge_Rec3_R[1:end-1]; Edge_Rec3_T[1:end-1]; Edge_Rec3_L[1:end]]
    Edge_Rec1           =   [Edge_Rec1_R[1:end-1]; Edge_Rec1_T[1:end-1]; Edge_Rec1_L[1:end]]
    A_Rec1_B, B_Rec1_B  =   loftlinear(Edge_Rec1, Edge_Rec3 ; num_steps=n_off, close_loop=false, face_type=:quad)

    # CUT REC - EDGE - LTCo
    Edge_Rec3                   =   [Edge_Rec3_B[1:end-1]; Edge_Rec3_R[1:end]]
    Edge_Rec1                   =   [Edge_Rec1_B[1:end-1]; Edge_Rec1_R[1:end]]
    A_Rec1_LTCo, B_Rec1_LTCo    =   loftlinear(Edge_Rec1, Edge_Rec3 ; num_steps=n_off, close_loop=false, face_type=:quad)

    # CUT REC - EDGE - LBCo
    Edge_Rec3                   =   [Edge_Rec3_R[1:end-1]; Edge_Rec3_T[1:end]]
    Edge_Rec1                   =   [Edge_Rec1_R[1:end-1]; Edge_Rec1_T[1:end]]
    A_Rec1_LBCo, B_Rec1_LBCo    =   loftlinear(Edge_Rec1, Edge_Rec3 ; num_steps=n_off, close_loop=false, face_type=:quad)

    # CUT REC - EDGE - RTCo
    Edge_Rec3                   =   [Edge_Rec3_L[1:end-1]; Edge_Rec3_B[1:end]]
    Edge_Rec1                   =   [Edge_Rec1_L[1:end-1]; Edge_Rec1_B[1:end]]
    A_Rec1_RTCo, B_Rec1_RTCo    =   loftlinear(Edge_Rec1, Edge_Rec3 ; num_steps=n_off, close_loop=false, face_type=:quad)

    # CUT REC - EDGE - RBCo
    Edge_Rec3                   =   [Edge_Rec3_T[1:end-1]; Edge_Rec3_L[1:end]]
    Edge_Rec1                   =   [Edge_Rec1_T[1:end-1]; Edge_Rec1_L[1:end]]
    A_Rec1_RBCo, B_Rec1_RBCo    =   loftlinear(Edge_Rec1, Edge_Rec3 ; num_steps=n_off, close_loop=false, face_type=:quad)

    # CUT REC - CHOOSE

    if     Edge == :I 
    A_Rec1, B_Rec1 = A_Rec1_I, B_Rec1_I
    elseif Edge == :R 
    A_Rec1, B_Rec1 = A_Rec1_R, B_Rec1_R
    elseif Edge == :L  
    A_Rec1, B_Rec1 = A_Rec1_L, B_Rec1_L
    elseif Edge == :T  
    A_Rec1, B_Rec1 = A_Rec1_T, B_Rec1_T
    elseif Edge == :B  
    A_Rec1, B_Rec1 = A_Rec1_B, B_Rec1_B
    elseif Edge == :LTCo
    A_Rec1, B_Rec1 = A_Rec1_LTCo, B_Rec1_LTCo
    elseif Edge == :LBCo
    A_Rec1, B_Rec1 = A_Rec1_LBCo, B_Rec1_LBCo
    elseif Edge == :RTCo
    A_Rec1, B_Rec1 = A_Rec1_RTCo, B_Rec1_RTCo
    elseif Edge == :RBCo
    A_Rec1, B_Rec1 = A_Rec1_RBCo, B_Rec1_RBCo
    else
        error("Unknown edge: $Edge")
    end


    # CUT REC - FINAL LOFT 
    B_combined      = [B_blank; B_Rec1]
    A_combined      = [A_blank;
                    [f .+ length(B_blank) for f in A_Rec1]]

    A_final, B_final= mergevertices(A_combined, B_combined)      

    return (A_final, B_final)
end

function remove_elements(B, A, points)
    node_ids = [
        argmin([sum(abs2, node - p) for node in B])
        for p in points
    ]

    node_ids = unique(node_ids)

    keep = [
        !any(id -> id in Tuple(f), node_ids)
        for f in A
    ]

    return B, A[keep], node_ids
end


function find_corner(edge; angle_tol_deg=30.0)
    corner_ids = Int[]

    for i in 2:length(edge)-1
        v1 = edge[i] - edge[i-1]
        v2 = edge[i+1] - edge[i]

        if norm(v1) < eps() || norm(v2) < eps()
            continue
        end

        c = dot(v1, v2) / (norm(v1) * norm(v2))
        c = clamp(c, -1.0, 1.0)

        angle = acosd(c)

        if angle > angle_tol_deg
            push!(corner_ids, i)
        end
    end

    return corner_ids, edge[corner_ids]
end

function create_fillet(
    A2,
    B2;
    r0 = 2.0,
    plane = :XY,
    loc = :LT,
    Point0 = Point{3, Float64}(-15.0, 7.5, 0.0),
    n = 3)

    # =============================================================================
    # SCRIPT START
    # =============================================================================

    # INPUTS - DERIVED
    if plane == :YZ || plane == :XY
        if         loc == :RB
            off_els1            =     4
            off_els2            =     4
            elseif loc == :RT
            off_els1            =     4
            off_els2            =     4
            elseif loc == :LB
            off_els1            =     4
            off_els2            =     4
            elseif loc == :LT
            off_els1            =     4
            off_els2            =     4
            else
            error("error - check loc. Use :RB, :RT, :LB, or :LT")
        end
        else
        error("unsupported plane: $plane")
    end

    # BOUNDARY EDGES
    B2_Edg_idx                  =       boundaryedges(A2)
    B2_Edg                      =       B2[edges2curve(B2_Edg_idx)]

    # LOCATE CORNERS - for making fillet - Point0, Point1 (away from the arc cor1), Point2 (away from the arc cor2)
    Point0_idx                  =       argmin([sum(abs2, node - Point0) for node in B2_Edg])
    DistShouldbe                =       2.75*r0

    candidate_idx               =       (Point0_idx + 1):length(B2_Edg)
    relative_idx                =       argmin(abs(norm(B2_Edg[i] - Point0) - DistShouldbe) for i in candidate_idx)
    Point1_idx                  =       first(candidate_idx) + relative_idx - 1
    Point1                      =       B2_Edg[Point1_idx]

    candidate_idx               =       (Point0_idx - 1):-1:1
    relative_idx                =       argmin(abs(norm(B2_Edg[i] - Point0) - DistShouldbe) for i in candidate_idx)
    Point2_idx                  =       first(candidate_idx) - relative_idx + 1
    Point2                      =       B2_Edg[Point2_idx]

    # LAYER 1 - REMOVE
    Edge0                       =       B2_Edg[Point2_idx : Point1_idx]
    B2, A2, rem_idx2            =       remove_elements(B2, A2, Edge0[2:end-1])
    B2_Edg_idx                  =       boundaryedges(A2)
    B2_Edg                      =       B2[edges2curve(B2_Edg_idx)]

    # LAYER 2 - REMOVE 
    Point1_2_idx                =       argmin([sum(abs2, node - Point1) for node in B2_Edg])
    Point2_2_idx                =       argmin([sum(abs2, node - Point2) for node in B2_Edg])
    Edge0                       =       B2_Edg[Point2_2_idx+1:Point1_2_idx-1]
    B2, A2, rem_idx2            =       remove_elements(B2, A2, Edge0[2:end-1])
    B2_Edg_idx                  =       boundaryedges(A2)
    B2_Edg                      =       B2[edges2curve(B2_Edg_idx)]

    # LAYER 3 - REMOVE
    Point1_3_idx                =       argmin([sum(abs2, node - Point1) for node in B2_Edg])
    Point2_3_idx                =       argmin([sum(abs2, node - Point2) for node in B2_Edg])
    Edge0                       =       B2_Edg[Point2_3_idx+2:Point1_3_idx-2]
    B2, A2, rem_idx2            =       remove_elements(B2, A2, Edge0[2:end-1])
    B2_Edg_idx                  =       boundaryedges(A2)
    B2_Edg                      =       B2[edges2curve(B2_Edg_idx)]

    # LAYER 4 - REMOVE
    Point1_4_idx                =       argmin([sum(abs2, node - Point1) for node in B2_Edg])
    Point2_4_idx                =       argmin([sum(abs2, node - Point2) for node in B2_Edg])
    Edge0                       =       B2_Edg[Point2_4_idx+3:Point1_4_idx-3]
    B2, A2, rem_idx2            =       remove_elements(B2, A2, Edge0[2:end-1])
    B2_Edg_idx                  =       boundaryedges(A2)
    B2_Edg                      =       B2[edges2curve(B2_Edg_idx)]
    A2, B2                      =       mergevertices(A2, B2)

    # REC 1 - INNER
    if          plane == :YZ
        if     loc == :RB
            Point1R1            =       Point{3,Float64}(Point0[1], Point0[2] + 0*r0 , Point0[3] + 1*r0)
            Point2R1            =       Point{3,Float64}(Point0[1], Point0[2] - 1*r0 , Point0[3] + 1*r0)
            Point3R1            =       Point{3,Float64}(Point0[1], Point0[2] - 1*r0 , Point0[3] - 1*r0)
            Point4R1            =       Point{3,Float64}(Point0[1], Point0[2] + 1*r0 , Point0[3] - 1*r0)
            Point5R1            =       Point{3,Float64}(Point0[1], Point0[2] + 1*r0 , Point0[3] + 0*r0)

        elseif loc == :RT
            Point1R1            =       Point{3,Float64}(Point0[1], Point0[2] + 0*r0 , Point0[3] - 1*r0)
            Point2R1            =       Point{3,Float64}(Point0[1], Point0[2] - 1*r0 , Point0[3] - 1*r0)
            Point3R1            =       Point{3,Float64}(Point0[1], Point0[2] - 1*r0 , Point0[3] + 1*r0)
            Point4R1            =       Point{3,Float64}(Point0[1], Point0[2] + 1*r0 , Point0[3] + 1*r0)
            Point5R1            =       Point{3,Float64}(Point0[1], Point0[2] + 1*r0 , Point0[3] - 0*r0)

        elseif loc == :LB
            Point1R1            =       Point{3,Float64}(Point0[1], Point0[2] - 0*r0 , Point0[3] + 1*r0)
            Point2R1            =       Point{3,Float64}(Point0[1], Point0[2] + 1*r0 , Point0[3] + 1*r0)
            Point3R1            =       Point{3,Float64}(Point0[1], Point0[2] + 1*r0 , Point0[3] - 1*r0)
            Point4R1            =       Point{3,Float64}(Point0[1], Point0[2] - 1*r0 , Point0[3] - 1*r0)
            Point5R1            =       Point{3,Float64}(Point0[1], Point0[2] - 1*r0 , Point0[3] + 0*r0)

        elseif loc == :LT
            Point1R1            =       Point{3,Float64}(Point0[1], Point0[2] - 0*r0 , Point0[3] - 1*r0)
            Point2R1            =       Point{3,Float64}(Point0[1], Point0[2] + 1*r0 , Point0[3] - 1*r0)
            Point3R1            =       Point{3,Float64}(Point0[1], Point0[2] + 1*r0 , Point0[3] + 1*r0)
            Point4R1            =       Point{3,Float64}(Point0[1], Point0[2] - 1*r0 , Point0[3] + 1*r0)
            Point5R1            =       Point{3,Float64}(Point0[1], Point0[2] - 1*r0 , Point0[3] - 0*r0)
        else
            error("error - check function")
        end
        elseif  plane == :XY
        if     loc == :RB
            Point1R1            =       Point{3,Float64}(Point0[1] + 0*r0 , Point0[2] + 1*r0, Point0[3])
            Point2R1            =       Point{3,Float64}(Point0[1] - 1*r0 , Point0[2] + 1*r0, Point0[3])
            Point3R1            =       Point{3,Float64}(Point0[1] - 1*r0 , Point0[2] - 1*r0, Point0[3])
            Point4R1            =       Point{3,Float64}(Point0[1] + 1*r0 , Point0[2] - 1*r0, Point0[3])
            Point5R1            =       Point{3,Float64}(Point0[1] + 1*r0 , Point0[2] + 0*r0, Point0[3])

        elseif loc == :RT
            Point1R1            =       Point{3,Float64}(Point0[1] + 0*r0 , Point0[2] - 1*r0, Point0[3])
            Point2R1            =       Point{3,Float64}(Point0[1] - 1*r0 , Point0[2] - 1*r0, Point0[3])
            Point3R1            =       Point{3,Float64}(Point0[1] - 1*r0 , Point0[2] + 1*r0, Point0[3])
            Point4R1            =       Point{3,Float64}(Point0[1] + 1*r0 , Point0[2] + 1*r0, Point0[3])
            Point5R1            =       Point{3,Float64}(Point0[1] + 1*r0 , Point0[2] - 0*r0, Point0[3])

        elseif loc == :LB
            Point1R1            =       Point{3,Float64}(Point0[1] - 0*r0 , Point0[2] + 1*r0, Point0[3])
            Point2R1            =       Point{3,Float64}(Point0[1] + 1*r0 , Point0[2] + 1*r0, Point0[3])
            Point3R1            =       Point{3,Float64}(Point0[1] + 1*r0 , Point0[2] - 1*r0, Point0[3])
            Point4R1            =       Point{3,Float64}(Point0[1] - 1*r0 , Point0[2] - 1*r0, Point0[3])
            Point5R1            =       Point{3,Float64}(Point0[1] - 1*r0 , Point0[2] + 0*r0, Point0[3])

        elseif loc == :LT
            Point1R1            =       Point{3,Float64}(Point0[1] - 0*r0 , Point0[2] - 1*r0, Point0[3])
            Point2R1            =       Point{3,Float64}(Point0[1] + 1*r0 , Point0[2] - 1*r0, Point0[3])
            Point3R1            =       Point{3,Float64}(Point0[1] + 1*r0 , Point0[2] + 1*r0, Point0[3])
            Point4R1            =       Point{3,Float64}(Point0[1] - 1*r0 , Point0[2] + 1*r0, Point0[3])
            Point5R1            =       Point{3,Float64}(Point0[1] - 1*r0 , Point0[2] - 0*r0, Point0[3])
        else
            error("error - check function")
        end
        else
        error("unsupported plane: $plane")
    end

    # REC 2 - OUTER
    Point1R2_idx                =       argmin([sum(abs2, node - Point1) for node in B2_Edg])
    Point5R2_idx                =       argmin([sum(abs2, node - Point2) for node in B2_Edg])
    if          loc == :RB
            Edge0               =       B2_Edg[Point5R2_idx:Point1R2_idx]
        elseif  loc == :RT
            Edge0               =       reverse(B2_Edg[Point5R2_idx:Point1R2_idx])
        elseif  loc == :LB
            Edge0               =       reverse(B2_Edg[Point5R2_idx:Point1R2_idx])
        elseif  loc == :LT
            Edge0               =       B2_Edg[Point5R2_idx:Point1R2_idx]
        else
            error("error - check loc")
    end

    Point1R2                    =       Edge0[1]
    Point1R2_idx                =       argmin([sum(abs2, node - Point1R2) for node in Edge0])
    Point5R2                    =       Edge0[end]
    Point5R2_idx                =       argmin([sum(abs2, node - Point5R2) for node in Edge0])
    Ed1Cor_idx, Ed1Cor          =       find_corner(Edge0)
    Point2R2                    =       Ed1Cor[1]
    Point2R2_idx                =       argmin([sum(abs2, node - Point2R2) for node in Edge0])
    Point3R2                    =       Ed1Cor[2]
    Point3R2_idx                =       argmin([sum(abs2, node - Point3R2) for node in Edge0])
    Point4R2                    =       Ed1Cor[3]
    Point4R2_idx                =       argmin([sum(abs2, node - Point4R2) for node in Edge0])

    # REC 2 - REC 1 - Creating the Partial Edges for Rec 2 Rec 1 and Loft Rec 2 with Rec 1
    Edge1R2                     =       Edge0[Point1R2_idx:Point2R2_idx]
    Edge1R1                     =       [(1.0-w)*Point1R1 + w*Point2R1 for w in range(0.0, 1.0, length=length(Edge1R2))]
    A3, B3                      =       loftlinear(Edge1R1 , Edge1R2 ; num_steps = n, close_loop = false, face_type = :quad)

    Edge2R2                     =       Edge0[Point2R2_idx:Point3R2_idx]
    Edge2R1                     =       [(1.0-w)*Point2R1 + w*Point3R1 for w in range(0.0, 1.0, length=length(Edge2R2))]
    A4, B4                      =       loftlinear(Edge2R1 , Edge2R2 ; num_steps = n, close_loop = false, face_type = :quad)

    Edge3R2                     =       Edge0[Point3R2_idx:Point4R2_idx]
    Edge3R1                     =       [(1.0-w)*Point3R1 + w*Point4R1 for w in range(0.0, 1.0, length=length(Edge3R2))]
    A5, B5                      =       loftlinear(Edge3R1 , Edge3R2 ; num_steps = n, close_loop = false, face_type = :quad)

    Edge4R2                     =       Edge0[Point4R2_idx:Point5R2_idx]
    Edge4R1                     =       [(1.0-w)*Point4R1 + w*Point5R1 for w in range(0.0, 1.0, length=length(Edge4R2))]
    A6, B6                      =       loftlinear(Edge4R1 , Edge4R2 ; num_steps = n, close_loop = false, face_type = :quad)

    # LOFT REC 1 - ARC EDGE
    EdgeR1                      =       [Edge2R1 ; Edge3R1[2:end]]
    if          plane == :YZ
        if     loc == :RB
            B_Arc               =       [Point{3,Float64}(Point1R1[1],Point1R1[2] + r0 - r0*cos(t), Point1R1[3] - r0*sin(t))
                                        for t in range(0.0, π/2, length=length(EdgeR1))]
        elseif loc == :RT
            B_Arc               =       [Point{3,Float64}(Point1R1[1],Point1R1[2] + r0 - r0*cos(t), Point1R1[3] - r0*sin(t))
                                        for t in range(2*π, 3*π/2, length=length(EdgeR1))]
        elseif loc == :LB
            B_Arc               =       [Point{3,Float64}(Point1R1[1],Point1R1[2] - r0 + r0*cos(t), Point1R1[3] - r0*sin(t))
                                        for t in range(0.0, π/2, length=length(EdgeR1))]

        elseif loc == :LT
            B_Arc               =       [Point{3,Float64}(Point1R1[1], Point1R1[2] - r0 + r0*cos(t), Point1R1[3] - r0*sin(t))
                                        for t in range(2π, 3π/2, length=length(EdgeR1))]
        else
            error("error - check function")
        end
        elseif  plane == :XY
            if     loc == :RB
            B_Arc               =       [Point{3,Float64}(Point1R1[1] + r0 - r0*cos(t), Point1R1[2] - r0*sin(t), Point1R1[3])
                                        for t in range(0.0, π/2, length=length(EdgeR1))]
        elseif loc == :RT
            B_Arc               =       [Point{3,Float64}(Point1R1[1] + r0 - r0*cos(t), Point1R1[2] - r0*sin(t), Point1R1[3])
                                        for t in range(2*π, 3*π/2, length=length(EdgeR1))]
        elseif loc == :LB
            B_Arc               =       [Point{3,Float64}(Point1R1[1] - r0 + r0*cos(t), Point1R1[2] - r0*sin(t), Point1R1[3])
                                        for t in range(0.0, π/2, length=length(EdgeR1))]
        elseif loc == :LT
            B_Arc               =       [Point{3,Float64}(Point1R1[1] - r0 + r0*cos(t), Point1R1[2] - r0*sin(t), Point1R1[3])
                                        for t in range(2π, 3π/2, length=length(EdgeR1))]
        else
            error("error - check function")
        end
        else
        error("unsupported plane: $plane")
    end

    A_Arc1, B_Arc1              =       loftlinear(B_Arc, EdgeR1; num_steps=length(Edge1R1), close_loop=false, face_type=:quad)

    # FINAL
    if     loc == :RB
        B_comb1                 =       [B3; B4; B5; B6; B_Arc1]
        A_comb1                 =       [A3;
                                        [(f) .+ (length(B3)) for f in A4];
                                        [(f) .+ (length(B3) + length(B4)) for f in A5];
                                        [(f) .+ (length(B3) + length(B4) + length(B5)) for f in A6];
                                        [(f) .+ (length(B3) + length(B4) + length(B5) + length(B6)) for f in A_Arc1]]
    elseif loc == :RT
        B_comb1                 =       [B3; B4; B5; B6; B_Arc1]
        A_comb1                 =       [[reverse(f) for f in A3];
                                        [reverse(f) .+ (length(B3)) for f in A4];
                                        [reverse(f) .+ (length(B3) + length(B4)) for f in A5];
                                        [reverse(f) .+ (length(B3) + length(B4) + length(B5)) for f in A6];
                                        [reverse(f) .+ (length(B3) + length(B4) + length(B5) + length(B6)) for f in A_Arc1]]
    elseif loc == :LB
        B_comb1                 =       [B3; B4; B5; B6; B_Arc1]
        A_comb1                 =       [[reverse(f) for f in A3];
                                        [reverse(f) .+ (length(B3)) for f in A4];
                                        [reverse(f) .+ (length(B3) + length(B4)) for f in A5];
                                        [reverse(f) .+ (length(B3) + length(B4) + length(B5)) for f in A6];
                                        [reverse(f) .+ (length(B3) + length(B4) + length(B5) + length(B6)) for f in A_Arc1]]

    elseif loc == :LT
        B_comb1                 =       [B3; B4; B5; B6; B_Arc1]
        A_comb1                 =       [A3;
                                        [(f) .+ (length(B3)) for f in A4];
                                        [(f) .+ (length(B3) + length(B4)) for f in A5];
                                        [(f) .+ (length(B3) + length(B4) + length(B5)) for f in A6];
                                        [(f) .+ (length(B3) + length(B4) + length(B5) + length(B6)) for f in A_Arc1]]
    else
        error("error - check function")
    end

    A_comb1, B_comb1        =   mergevertices(A_comb1, B_comb1)
    B_comb2                 =   [B2; B_comb1]
    A_comb2                 =   [A2;
                                [f .+ (length(B2)) for f in A_comb1]]
    A_final, B_final        =   mergevertices(A_comb2, B_comb2)

    # =============================================================================
    # SCRIPT - END
    # =============================================================================
    return (A_final, B_final)
end

function combine_faces_nodes(Bs, As, directions)
    sam(f) = f
    rev(f) = reverse(f)

    funcs = Dict(
        :sam => sam,
        :rev => rev,
    )

    offsets = cumsum([0; [length(B) for B in Bs[1:end-1]]])

    B = [
        b
        for B_part in Bs
        for b in B_part
    ]

    A = [
        funcs[directions[i]](f .+ offsets[i])
        for i in eachindex(As)
        for f in As[i]
    ]

    return A, B
end

function copy_nSections(A2, B2; nSections = 4, ShiftX = 200.0, ShiftY = 200.0, ShiftZ = 300.0,)
        A3                  =       Vector{QuadFace{Int64}}(undef,length(A2)*nSections)
        B3                  =       Vector{Point{3,Float64}}(undef,length(B2)*nSections)
        i_f = 1
        i_v = 1
        s = 0
        for q in 1:nSections
            A3[i_f:i_f+length(A2)-1]  = [f.+s for f in A2] 
            B3[i_v: i_v+length(B2)-1] = [Point{3, Float64}(v[1]- (q-1)*ShiftX, v[2] - (q-1)*ShiftY , v[3] - (q-1)*ShiftZ) for v in B2] 
            i_f += length(A2) 
            i_v += length(B2) 
            s   += length(B2)
        end
    return A3, B3
end # TO DO - for XY and YZ


function add_teeth(A1, B1; Teeth_wd = 20.0, Teeth_ht =20.0, Point0 = Point{3,Float64}(-25.0,  0.0, 0.0), Plane = :XY, Edge = :L,  off= 2.0,)
    ### SCRIPT - START
    if            Plane == :XY
        Point1_TnB              =    Point0 + [-(0.5*Teeth_wd), 0.0, 0.0]
        Point2_TnB              =    Point0 + [+(0.5*Teeth_wd), 0.0, 0.0]
        Point1off_TnB           =    Point1_TnB + [-off, 0.0, 0.0]
        Point2off_TnB           =    Point2_TnB + [+off, 0.0, 0.0]

        Point1_RnL              =    Point0 + [0.0, -(0.5*Teeth_wd), 0.0]
        Point2_RnL              =    Point0 + [0.0, +(0.5*Teeth_wd), 0.0]
        Point1off_RnL           =    Point1_RnL  + [0.0, -off, 0.0]
        Point2off_RnL           =    Point2_RnL  + [0.0, +off, 0.0]

    elseif        Plane == :YZ

        Point1_TnB              =    Point0 + [0.0, -(0.5*Teeth_wd), 0.0]
        Point2_TnB              =    Point0 + [0.0, +(0.5*Teeth_wd), 0.0]
        Point1off_TnB           =    Point1_TnB + [0.0, -off, 0.0]
        Point2off_TnB           =    Point2_TnB + [0.0, +off, 0.0]

        Point1_RnL              =    Point0 + [0.0, 0.0, -(0.5*Teeth_wd)]
        Point2_RnL              =    Point0 + [0.0, 0.0, +(0.5*Teeth_wd)]
        Point1off_RnL           =    Point1_RnL  + [0.0, 0.0, -off]
        Point2off_RnL           =    Point2_RnL  + [0.0, 0.0, +off]
        else
        error("error")
    end

    if            Edge ==  :T || Edge ==:B
        Point1      = Point1_TnB    
        Point2      = Point2_TnB
        Point1off   = Point1off_TnB
        Point2off   = Point2off_TnB
    elseif        Edge ==  :R || Edge ==:L
        Point1      = Point1_RnL    
        Point2      = Point2_RnL
        Point1off   = Point1off_RnL
        Point2off   = Point2off_RnL
        else
        error("error")
    end

    B1_Edg_idx          =   boundaryedges(A1)
    B1_Edg              =   B1[edges2curve(B1_Edg_idx)]
    Point1off_idx       =   argmin([sum(abs2, node - Point1off) for node in B1_Edg])
    Point1off           =   B1_Edg[Point1off_idx]
    Point2off_idx       =   argmin([sum(abs2, node - Point2off) for node in B1_Edg])
    Point2off           =   B1_Edg[Point2off_idx]

    # LAYER 1 - REMOVE
    i1, i2              =   minmax(Point1off_idx, Point2off_idx)
    EdgeL1              =   B1_Edg[i1 : i2]
    B1, A1, rem_idx1    =   remove_elements(B1, A1, EdgeL1[2:end-1])
    B1_Edg_idx          =   boundaryedges(A1)
    B1_Edg              =   B1[edges2curve(B1_Edg_idx)]
    A1, B1              =   mergevertices(A1, B1)
    EdgeL1              =   B1_Edg[i1 + 1 : i2 + 1]

    off_pts             =   2
    n_target            =   length(EdgeL1)
    n_side              =   off_pts
    n_middle            =   n_target - 2*n_side + 2

    @assert n_middle >= 2 "EdgeL1 does not have enough points"

    EdgeL2P1            =   [(1.0-w)*Point1off + w*Point1 for w in range(0.0, 1.0, length=n_side )]
    EdgeL2P2            =   [(1.0-w)*Point1    + w*Point2 for w in range(0.0, 1.0, length=n_middle)]
    EdgeL2P3            =   [(1.0-w)*Point2    + w*Point2off for w in range(0.0, 1.0, length=n_side )]
    EdgeL2              =   [EdgeL2P1[1:end-1]; EdgeL2P2[1:end-1] ;EdgeL2P3]
    if Plane == :XY
        EdgeL1_TnB = sort(EdgeL1; by = p -> p[1])
        EdgeL2_TnB = sort(EdgeL2; by = p -> p[1])
        EdgeL1_RnL = sort(EdgeL1; by = p -> p[2])
        EdgeL2_RnL = sort(EdgeL2; by = p -> p[2])

    elseif Plane == :YZ
        EdgeL1_TnB = sort(EdgeL1; by = p -> p[2])
        EdgeL2_TnB = sort(EdgeL2; by = p -> p[2])
        EdgeL1_RnL = sort(EdgeL1; by = p -> p[3])
        EdgeL2_RnL = sort(EdgeL2; by = p -> p[3])
        else
            error("error")
    end


    if          Edge ==:T || Edge ==:B
        EdgeL1 = EdgeL1_TnB
        EdgeL2 = EdgeL2_TnB
    elseif      Edge ==:R || Edge ==:L
        EdgeL1 = EdgeL1_RnL
        EdgeL2 = EdgeL2_RnL
        else
            error("error")
    end

    A2, B2              =    loftlinear(EdgeL1 , EdgeL2 ; num_steps = 2, close_loop = false, face_type = :quad)

    EdgeL2              =    EdgeL2[off_pts:end-off_pts+1]

    if Plane == :XY
        EdgeL3_T        =    [Point{3, Float64}(v[1], v[2] + Teeth_ht, v[3]) for v in EdgeL2]
        EdgeL3_B        =    [Point{3, Float64}(v[1], v[2] - Teeth_ht, v[3]) for v in EdgeL2]
        EdgeL3_R        =    [Point{3, Float64}(v[1] + Teeth_ht, v[2] , v[3]) for v in EdgeL2]
        EdgeL3_L        =    [Point{3, Float64}(v[1] - Teeth_ht, v[2] , v[3]) for v in EdgeL2]
    elseif Plane == :YZ
        EdgeL3_T        =    [Point{3, Float64}(v[1], v[2] , v[3]+ Teeth_ht) for v in EdgeL2]
        EdgeL3_B        =    [Point{3, Float64}(v[1], v[2] , v[3]- Teeth_ht) for v in EdgeL2]
        EdgeL3_R        =    [Point{3, Float64}(v[1] , v[2]+ Teeth_ht , v[3]) for v in EdgeL2]
        EdgeL3_L        =    [Point{3, Float64}(v[1] , v[2]- Teeth_ht , v[3]) for v in EdgeL2]
        else
         error("error")
    end

    if Edge == :T
        EdgeL3 = EdgeL3_T
    elseif Edge == :B
        EdgeL3 = EdgeL3_B
    elseif Edge == :R
        EdgeL3 = EdgeL3_R
    elseif Edge == :L
        EdgeL3 = EdgeL3_L
        else
        error("error")
    end

    A3, B3              =   loftlinear(EdgeL2 , EdgeL3 ; num_steps = 5, close_loop = false, face_type = :quad)

    A4_T, B4_T          = combine_faces_nodes([B1, B2, B3],[A1, A2, A3],[:sam, :rev, :rev])
    A4_B, B4_B          = combine_faces_nodes([B1, B2, B3],[A1, A2, A3],[:sam, :sam, :sam])
    A4_R, B4_R          = combine_faces_nodes([B1, B2, B3],[A1, A2, A3],[:sam, :sam, :sam])
    A4_L, B4_L          = combine_faces_nodes([B1, B2, B3],[A1, A2, A3],[:sam, :rev, :rev])
    if Edge == :T
        A4, B4 = A4_T, B4_T
    elseif Edge == :B
        A4, B4 = A4_B, B4_B
    elseif Edge == :R
        A4, B4 = A4_R, B4_R
    elseif Edge == :L
        A4, B4 = A4_L, B4_L
        else
        error("error")
    end

    # FINAL 
    A_final, B_final     =   mergevertices(A4, B4)

    return A_final, B_final
end
    
# TO DO - add hte partiition surface function, split it in XY and YZ position, and  set the example.