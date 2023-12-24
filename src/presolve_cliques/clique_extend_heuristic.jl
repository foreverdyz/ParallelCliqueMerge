#clique_extend_heuristic.jl

include("clique_extend.jl")

function clq_extend_parallel_heuristic(set::Vector{Vector{Int64}}, cg::SparseMatrixCSC, numthreads::Int64)
    m = length(set)
    (m < 1) && (return set, set)
    width = ceil(Int64, m / numthreads)
    con_set_g = [Int64[] for _ in 1:length(set)]
    lazy_set_g = [Vector{Int64}[] for _ in 1:length(set)]
    @threads for id in 1 : numthreads
        let
            for i in 1 : width
                (i % 2 == 1) ? (j = (i - 1) * numthreads + id) : (j = i * numthreads - id + 1)
                if j <= m
                    #C is cliques set, and pos is the index, where C[pos] is the longest extended cliques
                    C, pos = _extend_unique_clique(set[j], cg)
                    #send results back the main core
                    con_set_g[j] = C[pos]
                    #remove the pos term
                    local lazy_set = C[[k for k = 1:length(C) if k != pos]]
                    (length(lazy_set) > 0) && (lazy_set_g[j] = lazy_set)
                end
            end
        end
    end
    return con_set_g, lazy_set_g
end
