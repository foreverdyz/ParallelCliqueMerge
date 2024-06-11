#detect_cliques.jl

using BangBang
using Base.Threads

#detect cliques from conflicted knapsacks
function detect_cliques(knapsack_set::Vector{Tuple{Any, Any, Any}}, time_limit::Real)
    m = length(knapsack_set)
    clique_set = Vector{Int64}[]
    main_cliques = Vector{Int64}[]
    start_time = time();
    for j in 1:m
        if time() - start_time > time_limit
            break
        end
        local C = _unique_cliques_find(knapsack_set[j][1], knapsack_set[j][2], knapsack_set[j][3])
        if C != false
            push!(main_cliques, C[1])
            (length(C) > 1) && (append!!(clique_set, C[2:end]))
        end
    end
    return main_cliques, clique_set
end

function detect_cliques_parallel(knapsack_set::Vector{Tuple{Any, Any, Any}}, numthreads::Int64, time_limit::Real)
    m = length(knapsack_set)
    numthreads = min(m, numthreads)
    numthreads = max(1, numthreads)
    width = floor(Int64, m / numthreads)
    spl = SpinLock()
    clique_set_g = Vector{Int64}[]
    main_cliques_g = Vector{Int64}[]
    start_time = time();
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
                if time() - start_time > time_limit - 0.1 #(0.1 seconds for reduction of parallel computing)
                    break
                end
                local C = _unique_cliques_find(knapsack_set[j][1], knapsack_set[j][2], knapsack_set[j][3])
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

function _unique_cliques_find(a_index::Vector{Int64}, a_val::AbstractVector, b::Real)
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
        (a_val[a_sort[pos]] + a_val[a_sort[pos+1]] > b + 0.0001) ? (en = pos) : (st = pos)
    end
    #determine the pos
    (a_val[a_sort[st]] + a_val[a_sort[en]] > b + 0.0001) ? (pos = st) : (pos = en)
    (a_val[a_sort[pos]] + a_val[a_sort[pos+1]] <= b - 0.0001) && (println(false))

    #the main clique
    C = a_index[a_sort[pos:end]]

    (pos > 1) ? (index = pos - 1) : (return [C])
    #find additional maximum cliques
    C_full = [C]
    pos += 1
    while index > 0
        (pos > length(a_index)) && (return C_full)
        if a_val[a_sort[pos]] + a_val[a_sort[index]] > b + 0.0001
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