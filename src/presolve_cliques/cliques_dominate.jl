#cliques_dominate.jl

using SparseArrays

#check domination between cliques with 1 thread (serial)
function dominate_cliques(clique_set)
    #initialize an all-one list to denote whether one clique has been dominated or not
    cliques_list = zeros(length(clique_set))
    #the last clique does not need to check
    iter = 0
    for i in 1:(length(clique_set)-1)
        #only consider non-dominated cliques
        if cliques_list[i] < 1
            #local x = findnz(clique_set[i])[1]
            local x = clique_set[i]
            #check all non-checked cliques
            for j in (i+1):length(clique_set)
                #only check non-dominated cliques
                if cliques_list[j] < 1
                    #local y = findnz(clique_set[j])[1]
                    iter += 1
                    local y = clique_set[j]
                    if length(x) >= length(y)
                        #i dominates j
                        (issubset(y, x)) && (cliques_list[j] = 1)
                    else
                        #j dominates i
                        (issubset(x, y)) && (cliques_list[i] = 1)
                    end
                end
            end
        end
    end
    #output all non-dominated cliques
    return [clique_set[index] for index in 1:length(clique_set) if cliques_list[index] < 1]
end

#check domination between cliques with mulitple thread (parallel)
function dominate_cliques_parallel(clique_set, numthreads)
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
                            (issubset(y, x)) && (cliques_list[j] = 1)
                        else
                            #j dominates i
                            (issubset(x, y)) && (cliques_list[i] = 1)
                        end
                    end
                end
            end
        end
    end

    #output all non-dominated cliques
    return [clique_set[index] for index in 1:length(clique_set) if cliques_list[index] < 1]
end
