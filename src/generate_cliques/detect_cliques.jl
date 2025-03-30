#detect_cliques.jl

using BangBang, Base.Threads

"""
    detect_clqiues()
Detect clqiues from the set of knapsacks
"""
function detect_cliques(knapsack_set::Vector{Tuple{Any, Any, Any}}, numthreads::Int64, knapsack_size::Int64)
    m = length(knapsack_set)
    # Early exit if the input is empty
    if m < 1
        return Vector{Int64}[], Vector{Int64}[]
    end
    
    # Determine the number of threads to use
    numthreads = min(m, numthreads)
    
    # If only one thread, process sequentially
    if numthreads == 1
        main_cliques, clique_set = Vector{Int64}[], Vector{Int64}[]
        for j in 1:m
            if length(knapsack_set[j][1]) <= knapsack_size #avoid exceeding runtime
                C = _unique_cliques_find(knapsack_set[j][1], knapsack_set[j][2], knapsack_set[j][3])
                if C !== false
                    push!(main_cliques, C[1])
                    (length(C) > 1) && append!!(clique_set, C[2:end])
                end
	    end
        end
        return main_cliques, clique_set
    end
    
    # Parallel processing for more than one thread
    width = div(m, numthreads)
    spl = SpinLock()  # Ensure thread-safe access to shared variables
    main_cliques_g = [Int64[] for _ in 1:m]
    clique_set_g = Vector{Int64}[]

    @threads for id in 1:numthreads
        local_clique_set = Vector{Int64}[]  # Store local clique sets per thread

        # Determine the range of data this thread will process
        st = (id - 1) * width + 1
        en = (id == numthreads) ? m : min(id * width, m)

        # Process knapsack set in the thread's range
        for j in st:en
            if length(knapsack_set[j][1]) <= knapsack_size
                C = _unique_cliques_find(knapsack_set[j][1], knapsack_set[j][2], knapsack_set[j][3])
                if C !== false
                    main_cliques_g[j] = C[1]
                    if length(C) > 1
                        append!!(local_clique_set, C[2:end])
                    end
                end
            end
        end

        # Thread-safe merging of results
        lock(spl)
        append!!(clique_set_g, local_clique_set)
        unlock(spl)
    end

    # Remove empty clique entries
    main_cliques_g = filter!(x -> !isempty(x), main_cliques_g)
    return main_cliques_g, clique_set_g
end

"""
    _unique_clique_find()
Detect clqiues from one knapsack.
"""
function _unique_cliques_find(a_index::Vector{Int64}, a_val::AbstractVector, b::Real)
    # Sort `a_val` in non-decreasing order and store the indices
    a_sort = sortperm(a_val, rev=false)

    # Check boundary condition: if the sum of the two largest values is less than or equal to `b`
    if a_val[a_sort[end]] + a_val[a_sort[end-1]] <= b
        return false  # No clique exists for this constraint
    end

    #find clique C (the basic one)
    #find pivot by binary search
    st = 1
    en = length(a_index)
    while st + 1 < en
        pos = ceil(Int, (en + st) / 2)
        (a_val[a_sort[pos]] + a_val[a_sort[pos+1]] > b + 0.00001) ? (en = pos) : (st = pos)
    end
    #determine the pos
    pos = (a_val[a_sort[st]] + a_val[a_sort[en]] > b + 0.00001) ? st : en

    #the main clique
    C = a_index[a_sort[pos:end]]

    (pos > 1) ? (index = pos - 1) : (return [C])
    #find additional maximum cliques
    C_full = [C]
    pos += 1
    while index > 0
        (pos > length(a_index)) && (return C_full)
        # If the next pair exceeds `b`, create a new maximal clique
        if a_val[a_sort[pos]] + a_val[a_sort[index]] > b + 0.00001
            c = vcat(a_index[a_sort[pos:end]], a_index[a_sort[index]])
            push!(C_full, c)
            index -= 1
        else
            pos += 1
        end
    end
    return C_full #the first one is the first maximal clique
end
