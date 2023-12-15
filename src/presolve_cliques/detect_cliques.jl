#detect_cliques.jl

using SparseArrays

#find cliques from knapsack constraints
function detect_cliques(knapsack_set::Vector{Tuple{Any, Any}})
    m = length(knapsack_set)
    #initialize two vectors to store cliques detected from knapsack constraints
    #main_cliques includes cliques will be added to the original model
    #clique_set includes cliques will be used as lazy cuts (user cuts)
    clique_set = Vector{Int64}[]
    main_cliques = Vector{Int64}[]
    for j in 1:m
        local C = _unique_cliques_find(knapsack_set[j][1], knapsack_set[j][2])
        if C != false
            push!(main_cliques, C[1])
            (length(C) > 1) && (append!!(clique_set, C[2:end]))
        end
    end
    return main_cliques, clique_set
end

#for one knapsack constraint, find cliques
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

#also provide method to detect cliques in parallel 
function detect_cliques_parallel(knapsack_set::Vector{Tuple{Any, Any}}, numthreads::Int64)
    m = length(knapsack_set)
    numthreads = min(m, numthreads)
    numthreads = max(1, numthreads)
    width = floor(Int64, m / numthreads)
    spl = SpinLock()
    clique_set_g = Vector{Int64}[]
    main_cliques_g = Vector{Int64}[]
    @threads for id in 1:numthreads
        let
            local st = (id-1) * width + 1
            if id == numthreads
                local en = m
            else
                local en = id * width
            end
            
            local clique_set = Vector{Int64}[]
            local main_cliques = Vector{Int64}[]
            for j in st:en
                local C = _unique_cliques_find(knapsack_set[j][1], knapsack_set[j][2])
                if C != false
                    push!(main_cliques, C[1])
                    (length(C) > 1) && (append!!(clique_set, C[2:end]))
                end
            end
            lock(spl)
            append!!(main_cliques_g, main_cliques)
            append!!(clique_set_g, clique_set)
            unlock(spl)
        end
    end
    return main_cliques_g, clique_set_g
end