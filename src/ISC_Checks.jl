# ==============================================================================
# File        : Checks.jl
# Project     : ISC.jl
# Author      : Kamran Ali Bhatti
#
# Description :
# This file contains mesh validation functions used to verify the integrity
# of nodes, element connectivity and mesh quality before exporting the Abaqus
# input (.inp) file.
# ==============================================================================

function check1_quad_valid(A, B)
    # Check # 1 : Does every element in A point to valid nodes in B?  
    n_nodes = size(B, 1)
    n_elems = size(A, 1)

    for elem_id in 1:n_elems
        elem = A[elem_id]

        for node_id in elem
            if node_id < 1 || node_id > n_nodes
                println("Error")
                println("Element ", elem_id, " uses invalid node ID ", node_id)
                println("Valid node IDs are 1 to ", n_nodes)
                error("Invalid node ID found")
            end
        end
    end
    println("Check # 1, Ok, every elem in A point to valid node in B.")
    return true
end

function check2_quad_unique(A, B)
    # Check # 2, Does every quad face have 4 different node IDs?
    n_nodes = size(B, 1)
    n_elems = size(A, 1)
    for elem_id in 1:n_elems
        elem = A[elem_id]

        if length(unique(elem)) !=4
            println("Error")
            println("Element", elem_id, " does not have 4 unique nodes.")
            println("Element connectivity: ",elem)
            error("Duplicate node found inside one quad element")
        end
    end
    println("Check # 2, Ok, every element in A have 4 unique node IDs.")
    return true
end

function check3_quad_nonzero(A, B)
    # Check # 3, non zero area. 
    min_area = 1e-10

    for elem_id in 1:n_elems
        elem = A[elem_id]
        p1 = B[elem[1]]
        p2 = B[elem[2]]
        p3 = B[elem[3]]
        p4 = B[elem[4]]

        area1 = 0.5 * norm(cross(p2 - p1, p3 - p1))
        area2 = 0.5 * norm(cross(p3 - p1, p4 - p1))

        quad_area = area1 + area2

        if quad_area < min_area
            println("Error")
            println("Element ", element_id, " has near-zero area.")
            println("Element connectivity: ", element)
            println("Area: ", quad_area)
            error("Collapsed quad element found")
        end
    end
    println("Check # 3, Ok, all elements in A have nonzero area.")
    return true
end

function check4_quad_connected(A, B)
    # Check #  4: Is the mesh one connected component?
    node_to_elements = Dict{Int, Vector{Int}}()

    for element_id in 1:n_elems
        element = A[element_id]

        for node_id in element
            if !haskey(node_to_elements, node_id)
                node_to_elements[node_id] = Int[]
            end

            push!(node_to_elements[node_id], element_id)
        end
    end

    visited_elements = Set{Int}()
    elements_to_visit = [1]

    while !isempty(elements_to_visit)
        current_element = pop!(elements_to_visit)

        if current_element in visited_elements
            continue
        end

        push!(visited_elements, current_element)

        element = A[current_element]

        for node_id in element
            neighbor_elements = node_to_elements[node_id]

            for neighbor_element in neighbor_elements
                if !(neighbor_element in visited_elements)
                    push!(elements_to_visit, neighbor_element)
                end
            end
        end
    end

    if length(visited_elements) != n_elems
        println("Error")
        println("Mesh is disconnected.")
        println("Connected elements: ", length(visited_elements))
        println("Total elements: ", n_elems)
        error("Disconnected mesh found")
    end

    println("Check # 4, Ok, : all mesh in A is one connected component.")
    return true
end