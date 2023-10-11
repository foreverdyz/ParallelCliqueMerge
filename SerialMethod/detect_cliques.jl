#detect_cliques.jl

using SparseArrays

function detect_cliques(knapsack_set::Vector{Tuple{Any, Any}})
    m = length(knapsack_set)
    clique_set = Vector{Int64}[]
    main_cliques = Vector{Int64}[]
    for j in 1:m
        local C = _unique_cliques_find(knapsack_set[j][1], knapsack_set[j][2])
        if C != false
            push!(main_cliques, C[1])
            (length(C) > 1) && (append!!(clique_set, C[2:end]))
        end
    end
    #check identical elements
    #main_cliques = _check_identical_oneset(main_cliques)
    #clique_set = _check_identical_twosets(main_cliques, clique_set)
    return main_cliques, clique_set
end

function _unique_cliques_find(a::SparseVector, b::Real)
    a_index, a_val = findnz(a)
    #sort to non-decreasing order, and cache ordered index to a_sort
    a_sort = sortperm(a_val, rev=false)

    #check boundary condition to determine whether there is a clique
    (a_val[a_sort[end]] + a_val[a_sort[end-1]] <= b) && return false#does not have clique in this constraint

    #find clique C (the basic one)
    #find pivot by binary search
    st = 1
    en = length(a_index)
    while st + 1 < en
        pos = ceil(Int, (en + st) / 2)
        (a_val[a_sort[pos]] + a_val[a_sort[pos+1]] > b) ? (en = pos) : (st = pos)
    end
    #determine the pos
    (a_val[a_sort[st]] + a_val[a_sort[en]] > b) ? (pos = st) : (pos = en)

    #the main clique
    C = a_index[a_sort[pos:end]]

    (pos > 1) ? (index = pos - 1) : (return [C])
    #find additional maximum cliques
    C_full = [C]
    pos += 1
    while index > 0
        (pos > length(a_index)) && (return C_full)
        if a_val[a_sort[pos]] + a_val[a_sort[index]] > b
            c = a_index[a_sort[pos:end]]
            push!(c, a_index[a_sort[index]])
            push!(C_full, c)
            index -= 1
        else
            pos += 1
        end
    end
    return C_full #the first one is the first maximal clique
end

function _check_identical_oneset(main_cliques::Vector{Vector{Int64}})
    l = length(main_cliques)
    id_list = spzeros(l)
    for (index, x) in enumerate(main_cliques)
        if id_list[index] < 1
            for i in (index+1) : l
                if length(x) == length(main_cliques[i])
                    (issubset(x, main_cliques[i])) && (id_list[i] = 1)
                end
            end
        end
    end
    return [main_cliques[index] for index in 1:l if id_list[index] < 1]
end

function _check_identical_twosets(main_cliques::Vector{Vector{Int64}}, clique_set::Vector{Vector{Int64}})
    l = length(clique_set)
    id_list = spzeros(l)
    for (index, x) in enumerate(clique_set)
        if id_list[index] < 1
            for k in main_cliques
                if length(x) == length(k)
                    (issubset(x, k)) && (id_list[index] = 1)
                end
            end
            for i in (index+1) : (l-1)
                if length(x) == length(clique_set[i])
                    (issubset(x, clique_set[i])) && (id_list[i] = 1)
                end
            end
        end
    end
    return [clique_set[index] for index in 1:l if id_list[index] < 1]
end