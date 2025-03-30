#clique_domination.jl

using SparseArrays
using Base.Threads

function clique_domination(clique_set::Vector{Vector{Int64}}, numthreads::Int64)
    if numthreads > 1
        return _dominate_cliques_parallel(clique_set)
    else
        return _dominate_cliques(clique_set)
    end
end

function _dominate_cliques_parallel(clique_set::Vector{Vector{Int64}})
    m = length(clique_set)
    #initialize an all-one list to denote whether one clique has been dominated or not
    cliques_list = BitVector(undef, m)
    fill!(cliques_list, false)
    
    # Use parallel loop with @threads
    @threads for i in 1:(m-1)
        # Only check non-dominated cliques
        if !cliques_list[i]
            x = clique_set[i]
            # Check all remaining cliques
            for j in (i+1):m
                if !cliques_list[j]
                    y = clique_set[j]
                    if length(x) >= length(y)
                        # i dominates j
                        if cmp(x, y)
                            cliques_list[j] = true
                        end
                    else
                        # j dominates i
                        if cmp(y, x)
                            cliques_list[i] = true
                            break  # No need to check more for clique i
                        end
                    end
                end
            end
        end
    end
    
    # Return only non-dominated cliques
    return [clique_set[i] for i in 1:m if !cliques_list[i]]
end

function _dominate_cliques(clique_set::Vector{Vector{Int64}})
    m = length(clique_set)
    cliques_list = BitVector(undef, m)
    fill!(cliques_list, false)

    for i in 1:(m-1)
        if !cliques_list[i]
            x = clique_set[i]
            for j in (i+1):m
                if !cliques_list[j]
                    y = clique_set[j]
                    if length(x) >= length(y)
                        if cmp(x, y)
                            cliques_list[j] = true
                        end
                    else
                        if cmp(y, x)
                            cliques_list[i] = true
                            break  # No need to check more for clique i
                        end
                    end
                end
            end
        end
    end

    # Return only non-dominated cliques
    return [clique_set[i] for i in 1:m if !cliques_list[i]]
end

function cmp(x,y)
    for i in y
        if !(i in x)
            return false
        end
    end
    return true
end