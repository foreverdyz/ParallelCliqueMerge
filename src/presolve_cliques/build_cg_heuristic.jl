#build_cg_heuristic.jl

include("build_cg.jl")

using Base.Threads
using SparseArrays
using ThreadsX


function build_cg_distributed_heuristic(set_pack::Vector{Vector{Int64}}, binary_number::Int64, numthreads::Int64)
    cg_init = _init_cg(binary_number)
    m = length(set_pack)
    (m < 1) && (return cg_init)
    numthreads = min(m, numthreads)
    width = ceil(Int64, m / numthreads)
    cg_unit = spzeros(Int64, 2*binary_number, 2*binary_number)
    cg_g = [cg_unit for _ in 1:numthreads]
    #set_pack should have been sorted!!!!!!
    @threads for id in 1:numthreads
        let
            local a = Int64[]
            local b = Int64[]
            local s = 0
            for i in 1 : width
                (i % 2 == 1) ? (j = (i - 1) * numthreads + id) : (j = i * numthreads - id + 1)
                if j <= m
                    for k in 2:length(set_pack[j])
                        #index-1 to avoid diag term
                        for t in 1:(k-1)
                            if set_pack[j][k] > set_pack[j][t]
                                push!(a, set_pack[j][k])
                                push!(b, set_pack[j][t])
                            else
                                push!(a, set_pack[j][t])
                                push!(b, set_pack[j][k])
                            end
                            s += 1
                        end
                    end
                end
            end
            push!(a, 2*binary_number)
            push!(b, 2*binary_number)
            c = ones(Int64, s)
            push!(c, 0)
            cg_g[id] = dropzeros!(sparse(a, b, c))
        end
    end
    (numthreads == 1) && (return cg_g[1] .| cg_init)
    L = ceil(Int64, log2(numthreads))
    K = numthreads
    for _ in 1:L
        local k = ceil(Int64, K/2)
        @threads for i in 1:k
            (i + k <= K) && (cg_g[i] = cg_g[i] .| cg_g[i+k])
        end
        K = k
    end
    return cg_g[1] .| cg_init
end

function update_cg_distributed_heuristic(set::Vector{Vector{Int64}}, binary_number::Int64, cg_init::SparseMatrixCSC, numthreads::Int64)
    m = length(set)
    (m < 1) && (return cg_init)
    numthreads = min(m, numthreads)
    width = ceil(Int64, m / numthreads)
    cg_unit = spzeros(Int64, 2*binary_number, 2*binary_number)
    cg_g = [cg_unit for _ in 1:numthreads]
    @threads for id in 1:numthreads
        let
            local a = Int64[]
            local b = Int64[]
            local s = 0
            for i in 1 : width
                (i % 2 == 1) ? (j = (i - 1) * numthreads + id) : (j = i * numthreads - id + 1)
                if j <= m
                    for k in 2:length(set[j])
                        #index-1 to avoid diag term
                        for t in 1:(k-1)
                            if set[j][k] > set[j][t]
                                push!(a, set[j][k])
                                push!(b, set[j][t])
                            else
                                push!(a, set[j][t])
                                push!(b, set[j][k])
                            end
                            s += 1
                        end
                    end
                end
            end
            push!(a, 2*binary_number)
            push!(b, 2*binary_number)
            c = ones(Int64, s)
            push!(c, 0)
            cg_g[id] = dropzeros!(sparse(a, b, c))
        end
    end
    (numthreads == 1) && (return cg_g[1] .| cg_init)
    L = ceil(Int64, log2(numthreads))
    K = numthreads
    for _ in 1:L
        local k = ceil(Int64, K/2)
        @threads for i in 1:k
            (i + k <= K) && (cg_g[i] = cg_g[i] .| cg_g[i+k])
        end
        K = k
    end
    return cg_g[1] .| cg_init
end

function heuristic_sort(set::Vector{Vector{Int64}})
    (length(set) < 2) && (return set) 
    set_length = [length(x) for x in set];
    set_ord = sortperm(set_length; alg = ThreadsX.MergeSort, rev = true);
    return [set[set_ord[i]] for i in 1:length(set)]
end
