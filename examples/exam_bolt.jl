using ISC
using ISC.Comodo
using ISC.Comodo.GLMakie

# INPUTS
# Center         =           (xc1, yc1, zc1) =     Point{3, Float64}(0.0 , 0.0 , 0.0)
# r2              =           6.0
# Th_BHead        =           5.3
# L_Thread        =           3.0+6.4+3.0          

A1, B1, E1 = comp_bolt(Center  = Point{3, Float64}(0.0 , 0.0 , 0.0), r_BHead = 6.0, Th_BHead = 5.3, r_Thread = 4.0, L_Thread = 12.4, n = 6)
A2, B2, E2 = comp_bolt(Center  = Point{3, Float64}(0.0 , 0.0 , 0.0), r_BHead = 6.0, Th_BHead = 5.3, r_Thread = 4.0, L_Thread = 18.0, n = 6)
A3, B3, E3 = comp_bolt(Center  = Point{3, Float64}(0.0 , 0.0 , 0.0), r_BHead = 6.0, Th_BHead = 5.3, r_Thread = 4.0, L_Thread = 30.0, n = 6)

### VISUALISATION
markersize = 18
linewidth  = 1.5
fig = Figure(size = (1800, 850))
# ax1 = Axis(fig[1, 1], xlabel = "X", ylabel = "Y")
ax1 = AxisGeom(fig[1, 1], title="Bolt1")
ax2 = AxisGeom(fig[1, 2], title="Bolt2")
ax3 = AxisGeom(fig[1, 3], title="Bolt3")

# scatter!(ax1, B_Bolt  , markersize = markersize *0.50    , color = :red)

meshplot!(ax1, A1, B1, strokewidth = 0.6, color = (:grey70, 0.75), strokecolor = (:black, 0.60), transparency = true)
meshplot!(ax2, A2, B2, strokewidth = 0.6, color = (:grey70, 0.75), strokecolor = (:black, 0.60), transparency = true)
meshplot!(ax3, A3, B3, strokewidth = 0.6, color = (:grey70, 0.75), strokecolor = (:black, 0.60), transparency = true)
# normalplot(ax1,A_Bolt,B_Bolt)

origin_line_length1 = 10
origin_line_length2 = 10
origin_line_length3 = 10

for ax in (ax1, ax2, ax3), (x, y, z, c) in (([-origin_line_length1, origin_line_length1], [0, 0], [0, 0], :red), ([0, 0], [-origin_line_length2, origin_line_length2], [0, 0], :green), ([0, 0], [0, 0], [-origin_line_length3, origin_line_length3], :blue))
    lines!(ax, x, y, z, color = c, linewidth = 1.5, linestyle = :dot, overdraw = true)
end

# wireframe!(ax1, GeometryBasics.Mesh(Vh, Eh), linewidth = 3, color = :red)
# normalplot(ax1,A1,B1)

fig

# out_dir = raw"C:\Users\23116524\OneDrive - National University of Ireland, Galway\1-PhD NUI Galway\AISC IDEAS AWARD\Pictures\03_Digital Platform"
# mkpath(out_dir)

# png_file = joinpath(out_dir, "exam_bolt.png")
# save(png_file, fig, px_per_unit = 8)

# println("Saved visualisation to: ", png_file)