#clique_extension.jl

using Base.Threads
using SparseArrays

# Main function for clique extension
function clique_extension(
        set::Vector{Vector{Int64}}, cg::SparseMatrixCSC, cg_trans::SparseMatrixCSC,
        numthreads::Int64,  clique_size::Int64, process_limit::Int64, time_limit::Real)
    if numthreads > 1
        return _clq_extend_parallel(set, cg, cg_trans, numthreads, clique_size, process_limit, time_limit)
    else
        return _clq_extend(set, cg, cg_trans, clique_size, process_limit, time_limit)
    end
end

function _clq_extend_parallel(
        set::Vector{Vector{Int64}}, cg::SparseMatrixCSC,  cg_trans::SparseMatrixCSC,
        numthreads::Int64, clique_size::Int64, process_limit_g::Int64, time_limit::Real)
    #set = shuffle(set_pack)
    #here we should shuffle the index set for set instead of the set itself to get a fixed results
    process_limit = ceil(Int64, process_limit_g/log2(numthreads+1));
    m = length(set)
    (m < 1) && (return set, [set])
    width = floor(Int64, m / numthreads)
    con_set_g = [Int64[] for _ in 1:length(set)]
    lazy_set_g = [Vector{Int64}[] for _ in 1:length(set)]
    index_set = shuffle(1:m)
    start_time = time();
    @threads for id in 1:numthreads
        let
            S = 0
            st = (id-1) * width + 1
            en = (id == numthreads) ? m : id * width
            
            for j in st:en
                i = index_set[j]
                if S > process_limit || length(set[i]) > clique_size || time() - start_time > time_limit
                  con_set_g[i] = set[i]
                else
                    C, pos, s = _extend_unique_clique(set[i], cg, cg_trans)
                    con_set_g[i] = C[pos]
                    local lazy_set = Vector{Int64}[]
                    for k in 1:length(C)
                        (k != pos) && (push!(lazy_set, C[k]))
                    end
                    (length(lazy_set) > 0) && (lazy_set_g[i] = lazy_set) 
                    S += s
                end
            end
        end
    end
    return con_set_g, lazy_set_g
end

function _clq_extend(set::Vector{Vector{Int64}}, cg::SparseMatrixCSC,  cg_trans::SparseMatrixCSC,
        clique_size::Int64, process_limit::Int64, time_limit::Real)
    m = length(set)
    con_set = Vector{Int64}[]
    lazy_set = Vector{Int64}[]
    S = 0
    start_time = time();
    for i in 1:m
        if S > process_limit || length(set[i]) > clique_size || time() - start_time > time_limit
            push!(con_set, set[i])
        else
            C, pos, s = _extend_unique_clique(set[i], cg, cg_trans)
            push!(con_set, C[pos])
            for i in 1:length(C)
                (i != pos) && (push!(lazy_set, C[i]))
            end
            S += s
        end
    end
    return con_set, [lazy_set]
end

function _extend_unique_clique(c::AbstractVector, cg::SparseMatrixCSC, cg_trans::SparseMatrixCSC,)
    #degree claculate
    degree_list = [sum(cg_trans[:, index]) + sum(cg[:, index]) for index in c]
    
    #s denote the total computation costs
    s = length(c)
    
    #get the smallest degree variables
    d = c[findmin(degree_list)[2]]
    #initialize an empty list for candidates
    candidate_list = Int64[]

    #find all potential candidates by for loop
    #for all conflict variables to d, there are two part as our represent way of CG
    #first part
    for index in findnz(cg_trans[:, d])[1]
        s += length(c) #costs to check index in c
        if !(index in c)
            local candidate_id = true
            #check whether variable index has been connected to all variables from c in CG
            for i in c
                s += 1
                if cg[index, i] + cg[i, index] < 1 
                    candidate_id = false
                    break
                end
            end
            #if yes, add variable index to the candidate list
            (candidate_id) && (push!(candidate_list, index))
        end
    end
    #second part
    for index in findnz(cg[:, d])[1]
        s += length(c)
        if !(index in c)
            local candidate_id = true
            for i in c
                s += 1
                if cg[index, i] + cg[i, index] < 1
                    candidate_id = false
                    break
                end
            end
            (candidate_id) && (push!(candidate_list, index))
        end
    end
    #if there is no candidate, just return the original clique
    (length(candidate_list) > 0) ? (extend_list = [[candidate_list[1]]]) : (return [c], 1, s)
    (s < 1_000_000) ? (extend_list = [[candidate_list[1]]]) : (return [c], 1, s)
    #now try to combine candidates to get stronger clique
    for i in candidate_list[2:end]
        #to check wether variable i cannot be combined to any other candidate
        #if yes, we need to generate a clique as original clique + variable i
        local isunique = true
        for (index, k) in enumerate(extend_list)
            #check whether i can be connected to the current clique
            local inthislist = true 
            for j in k
                s += 1
                if cg[i, j] + cg[j, i] < 1
                    inthislist = false
                    break
                end
            end
            if inthislist
                isunique = false
                push!(extend_list[index], i)
            end
        end
        (isunique) && (push!(extend_list, [i]))
        if s > 1_000_000
            break
        end
    end
    
    #find largest length term in extend_list
    pos = 0
    len = 0
    for (index, k) in enumerate(extend_list)
        (length(k) > len) && ((pos, len) = (index, length(k)))
    end
    
    return [append!!(k, c) for k in extend_list], pos, s
end
