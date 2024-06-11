#cg_construct.jl

using Random
using SparseArrays
using Base.Threads

"""
length_limit the size limit of clq we process; 
term_limit is the number of nonzero elements we process at most in one thread!!!
"""
function cg_construct(
        set_pack::Vector{Vector{Int64}}, set_pack_new::Vector{Vector{Int64}}, 
        main_cliques::Vector{Vector{Int64}}, clique_set::Vector{Vector{Int64}},
        org_to_bin::Dict{Int64, Int64}, numthreads::Int64,
        length_limit::Int64, term_limit::Int64
    )
    if numthreads > 1
        cg = _build_cg_distributed(set_pack, length(org_to_bin), numthreads, length_limit, term_limit);
        cg = _update_cg_distributed(set_pack_new, length(org_to_bin), cg, numthreads, length_limit, term_limit);
        cg = _update_cg_distributed(main_cliques, length(org_to_bin), cg, numthreads, length_limit, term_limit);
        cg = _update_cg_distributed(clique_set, length(org_to_bin), cg, numthreads, length_limit, term_limit);
        return min.(cg, 1)
    else
        cg = _build_cg(set_pack, length(org_to_bin), length_limit, term_limit);
        cg = _update_cg(set_pack_new, length(org_to_bin), cg, length_limit, term_limit);
        cg = _update_cg(main_cliques, length(org_to_bin), cg, length_limit, term_limit);
        cg = _update_cg(clique_set, length(org_to_bin), cg, length_limit, term_limit);
        return min.(cg, 1)
    end
end


function _build_cg_distributed(
    set::Vector{Vector{Int64}}, binary_number::Int64, numthreads::Int64,
    length_limit::Int64, term_limit::Int64
    )
    set_pack = shuffle(set)
    cg_init = _init_cg(binary_number)
    m = length(set_pack)
    (m < 1) && (return cg_init)
    numthreads = min(m, numthreads)
    numthreads = max(1, numthreads)
    width = floor(Int64, m / numthreads)
    spl = SpinLock()
    cg_unit = spzeros(Int64, 2*binary_number, 2*binary_number)
    cg_g = [cg_unit for _ in 1:numthreads]
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
            #here we want diag element keeps zero
                en_inner = Int64(min(length_limit, length(set_pack[j])))
                for i in 2:en_inner
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
                end
                if s > term_limit
                    break
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
    cg = cg_g[1] .| cg_init
    return cg
end

function _update_cg_distributed(
    set_pack::Vector{Vector{Int64}}, binary_number::Int64, 
    cg_init::SparseMatrixCSC, numthreads::Int64,
    length_limit::Int64, term_limit::Int64
    )
    set = shuffle(set_pack)
    m = length(set)
    (m < 1) && (return cg_init)
    numthreads = min(m, numthreads)
    numthreads = max(1, numthreads)
    width = floor(Int64, m / numthreads)
    spl = SpinLock()
    cg_unit = spzeros(Int64, 2*binary_number, 2*binary_number)
    cg_g = [cg_unit for _ in 1:numthreads]
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
            #here we want diag element keeps zero
                en_inner = Int64(min(length_limit, length(set[j])))
                for i in 2:en_inner
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
                end
                if s > term_limit
                    break
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
    cg = cg_g[1] .| cg_init
    return cg
end

function _build_cg(
    set::Vector{Vector{Int64}}, binary_number::Int64,
    length_limit::Int64, term_limit::Int64
    )
    set_pack = copy(set) #avoid changing set
    cg_init = _init_cg(binary_number)
    m = length(set_pack)
    (m < 1) && (return cg_init)
    cg = spzeros(2*binary_number, 2*binary_number)
    local a = Int64[]
    local b = Int64[]
    local s = 0
    for j in 1:m
        #here we want diag element keeps zero
        en_inner = Int64(min(length_limit, length(set_pack[j])))
        for i in 2:en_inner
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
        end
        if s > term_limit
            break
        end
    end
    push!(a, 2*binary_number)
    push!(b, 2*binary_number)
    c = ones(Int64, s)
    push!(c, 0)
    cg = dropzeros!(sparse(a, b, c)) 
    return cg .| cg_init
end

function _update_cg(
    set_pack::Vector{Vector{Int64}}, binary_number::Int64, cg_init::SparseMatrixCSC,
    length_limit::Int64, term_limit::Int64
    )
    set = copy(set_pack) #avoid changing set_pack
    m = length(set)
    (m < 1) && (return cg_init)
    cg = spzeros(2*binary_number, 2*binary_number)
    local a = Int64[]
    local b = Int64[]
    local s = 0
    for j in 1:m
    #here we want diag element keeps zero
        en_inner = Int64(min(length_limit, length(set[j])))
        for i in 2:en_inner
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
        end
        if s > term_limit
            break
        end
    end
    push!(a, 2*binary_number)
    push!(b, 2*binary_number)
    c = ones(Int64, s)
    push!(c, 0)
    cg = dropzeros!(sparse(a, b, c)) 
    return cg .| cg_init
end

function _init_cg(binary_number::Int64)
    #initialize an empty sparse matrix
    cg = spzeros(Int64, 2*binary_number, 2*binary_number)
    #x + bar(x) = 1, so x + bar(x) <= 1
    for i in 1:binary_number
        cg[i + binary_number, i] = 1
    end
    return cg
end