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
    cliques_list = zeros(m)
    #the last clique does not need to check
    @threads for i in 1:(m-1)
        #only consider non-dominated cliques
        if cliques_list[i] < 1
            let
                local x = clique_set[i]
                #check all non-checked cliques
                for j in (i+1):m
                    #only check non-dominated cliques
                    if cliques_list[j] < 1
                        local y = clique_set[j]
                        if length(x) >= length(y)
                            #i dominates j
                            (cmp(x, y)) && (cliques_list[j] = 1)
                        else
                            #j dominates i
                            (cmp(y, x)) && (cliques_list[i] = 1)
                        end
                    end
                end
            end
        end
    end
    
    #output all non-dominated cliques
    return [clique_set[index] for index in 1:length(clique_set) if cliques_list[index] < 1]
end

function _dominate_cliques(clique_set::Vector{Vector{Int64}})
    m = length(clique_set)
    #initialize an all-one list to denote whether one clique has been dominated or not
    cliques_list = zeros(m)
    #the last clique does not need to check
    for i in 1:(m-1)
        #only consider non-dominated cliques
        if cliques_list[i] < 1
            let
                local x = clique_set[i]
                #check all non-checked cliques
                for j in (i+1):m
                    #only check non-dominated cliques
                    if cliques_list[j] < 1
                        local y = clique_set[j]
                        if length(x) >= length(y)
                            #i dominates j
                            (cmp(x, y)) && (cliques_list[j] = 1)
                        else
                            #j dominates i
                            (cmp(y, x)) && (cliques_list[i] = 1)
                        end
                    end
                end
            end
        end
    end
    
    #output all non-dominated cliques
    return [clique_set[index] for index in 1:length(clique_set) if cliques_list[index] < 1]
end

function cmp(x,y)
    for i in y
        if !(i in x)
            return false
        end
    end
    return true
end