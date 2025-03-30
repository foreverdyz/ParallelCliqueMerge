#cg_construct.jl

using Random
using SparseArrays
using Base.Threads

"""
    cg_construct()
Construct conflict graph (CG) based on provided clqiues.
"""
function cg_construct(
        set_pack::Vector{Vector{Int64}}, set_pack_new::Vector{Vector{Int64}}, 
        main_cliques::Vector{Vector{Int64}}, clique_set::Vector{Vector{Int64}},
        org_to_bin::Dict{Int64, Int64}, numthreads::Int64,
        clique_size::Int64, process_limit::Int64, time_limit::Real
    )
    time_limit_here = time_limit/log2(numthreads+1);
    process_limit_here = ceil(Int64, process_limit/log2(numthreads+1));
    t = time();
    cg = _build_cg_distributed(set_pack, length(org_to_bin), numthreads, clique_size, process_limit_here, time_limit_here);
    t_used = time() - t;
    if time_limit_here - t_used > 0
        cg = _update_cg_distributed(set_pack_new, length(org_to_bin), cg, numthreads, clique_size, process_limit_here, time_limit_here - t_used);
    end
    t_used = time() - t;
    if time_limit_here - t_used > 0
        cg = _update_cg_distributed(main_cliques, length(org_to_bin), cg, numthreads, clique_size, process_limit_here, time_limit_here - t_used);
    end
    t_used = time() - t;
    if time_limit_here - t_used > 0
        cg = _update_cg_distributed(clique_set, length(org_to_bin), cg, numthreads, clique_size, process_limit_here, time_limit_here - t_used);
    end
    
    cg = min.(cg, 1)
    I, J, K = findnz(cg);
    cg_trans = sparse(J, I, K, 2*length(org_to_bin), 2*length(org_to_bin))
    return cg, cg_trans
end


function _build_cg_distributed(
        set::Vector{Vector{Int64}}, binary_number::Int64, numthreads::Int64,
        clique_size::Int64, process_limit::Int64, time_limit::Real
    )
    set_pack = shuffle(set)
    cg_init = _init_cg(binary_number)
    m = length(set_pack)
    (m < 1) && (return cg_init)
    numthreads = min(m, numthreads)
    width = div(m, numthreads)
    spl = SpinLock()
    cg_unit = spzeros(Int64, 2*binary_number, 2*binary_number)
    cg_g = [cg_unit for _ in 1:numthreads]
    start_time = time();
    @threads for id in 1:numthreads
        let
            local st = (id-1) * width + 1
            if id == numthreads
                local en = m
            else
                local en = id * width
            end
            local a = Int64[]
            local b = Int64[]
            local s = 0
            for j in st:en
                st_inner, en_inner = 1, length(set_pack[j])
                if en_inner > clique_size
                    st_inner = rand(1:(en_inner - clique_size + 1))
                    en_inner = st_inner + clique_size - 1
                end
                #here we want diag element keeps zero
                for i in st_inner+1 : en_inner
                    #index-1 to avoid diag term
                    for k in 1:(i-1)
                        if set_pack[j][i] > set_pack[j][k]
                            push!(a, set_pack[j][i])
                            push!(b, set_pack[j][k])
                        else
                            push!(a, set_pack[j][k])
                            push!(b, set_pack[j][i])
                        end
                        s += 1
                    end
                    if s > process_limit/numthreads || time() - start_time > time_limit #may change to others with larger memory like (log2(numthreads) + 1)
                        break
                    end
                end
            end
            c = ones(Int64, s)
            cg_g[id] = sparse(a, b, c, 2*binary_number, 2*binary_number) 
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
    cg = cg_g[1] .| cg_init
    return cg
end

function _update_cg_distributed(
        set_pack::Vector{Vector{Int64}}, binary_number::Int64, cg_init::SparseMatrixCSC, 
        numthreads::Int64, clique_size::Int64, process_limit::Int64, time_limit::Real
    )
    process_limit_here = process_limit - SparseArrays.nnz(cg_init)
    set = shuffle(set_pack)
    m = length(set)
    (m < 1) && (return cg_init)
    numthreads = min(m, numthreads)
    width = div(m, numthreads)
    spl = SpinLock()
    cg_unit = spzeros(Int64, 2*binary_number, 2*binary_number)
    cg_g = [cg_unit for _ in 1:numthreads]
    start_time = time();
    @threads for id in 1:numthreads
        let
            local st = (id-1) * width + 1
            if id == numthreads
                local en = m
            else
                local en = id * width
            end
            local a = Int64[]
            local b = Int64[]
            local s = 0
            for j in st:en
                st_inner, en_inner = 1, length(set[j])
                if en_inner > clique_size
                    st_inner = rand(1:(en_inner - clique_size + 1))
                    en_inner = st_inner + clique_size - 1
                end
                #here we want diag element keeps zero
                for i in st_inner+1 : en_inner
                    #index-1 to avoid diag term
                    for k in 1:(i-1)
                        if set[j][i] > set[j][k]
                            push!(a, set[j][i])
                            push!(b, set[j][k])
                        else
                            push!(a, set[j][k])
                            push!(b, set[j][i])
                        end
                        s += 1
                    end
                    if s > process_limit_here/numthreads || time() - start_time > time_limit #may change to others with larger memory like (log2(numthreads) + 1)
                        break
                    end
                end
            end
            c = ones(Int64, s)
            cg_g[id] = sparse(a, b, c, 2*binary_number, 2*binary_number)
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
    cg = cg_g[1] .| cg_init
    return cg
end

function _init_cg(binary_number::Int64)
    # Initialize an empty sparse matrix for the conflict graph
    cg = spzeros(Int64, 2*binary_number, 2*binary_number)
    # Add complementarity constraints (x + bar(x) <= 1)
    for i in 1:binary_number
        cg[i + binary_number, i] = 1
    end
    return cg
end
