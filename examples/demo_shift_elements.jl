using ISC
using GeometryBasics
# PART 1
V1 = [Point(0.0, 0.0, 0.0),Point(1.0, 0.0, 0.0),Point(1.0, 1.0, 0.0),Point(0.0, 1.0, 0.0),
      Point(0.0, 0.0, 1.0),Point(1.0, 0.0, 1.0), Point(1.0, 1.0, 1.0)Point(0.0, 1.0, 1.0)]

E1 = [(1,2,3,4,5,6,7,8)]
# PART 2
V2 = [Point(2.0, 0.0, 0.0),Point(3.0, 0.0, 0.0),Point(3.0, 1.0, 0.0),Point(2.0, 1.0, 0.0),
      Point(2.0, 0.0, 1.0),Point(3.0, 0.0, 1.0),Point(3.0, 1.0, 1.0),Point(2.0, 1.0, 1.0)]

E2 = [(1,2,3,4,5,6,7,8)]
# PART 3
V3 = [Point(4.0, 0.0, 0.0),Point(5.0, 0.0, 0.0),Point(5.0, 1.0, 0.0),Point(4.0, 1.0, 0.0),
     Point(4.0, 0.0, 1.0),Point(5.0, 0.0, 1.0),Point(5.0, 1.0, 1.0),Point(4.0, 1.0, 1.0)]
E3 = [(1,2,3,4,5,6,7,8)]
# ------------------------

# function shift_elements(E, offset)
#     [typeof(e)(ntuple(i -> Tuple(e)[i] + offset, length(Tuple(e)))) for e in E]
# end
# ------------------------
V = vcat(V1, V2, V3)
E = vcat(E1,
         shift_elements(E2, length(V1)),
         shift_elements(E3, length(V1) + length(V2)))
# ------------------------

# CHECK
# -------------------------------
println("Total Vertices = ", length(V))
println("Total Elements = ", length(E))
println("Nodes in first element = ", length(Tuple(E[1])))
for (i, e) in enumerate(E)
    println("Element $i => ", Tuple(e))
end