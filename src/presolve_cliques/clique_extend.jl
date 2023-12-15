#cliques_extend.jl

using BangBang
using Base.Threads
using SparseArrays


function clq_extend_parallel_new(set::Vector{Vector{Int64}}, cg::SparseMatrixCSC, numthreads::Int64)
    m = length(set)
    #width is the lower bound of the number of constraints for each thread
    width = floor(Int64, m / numthreads)
    lazy_set_g = [Vector{Int64}[] for _ in 1:length(set)]
    #parallel computation
    @threads for id in 1:numthreads
        let
            #determine indeces of constraints for each thread
            local st = (id - 1) * width + 1
            if id == numthreads
                local en = m
            else
                local en = id * width
            end

            #extend cliques
            for i in st:en
                #C is cliques set, and pos is the index, where C[pos] is the longest extended cliques
                C, pos = _extend_unique_clique(set[i], cg)
                #send results back the main core
                con_set_g[i] = C[pos]
                #remove the pos term
                local lazy_set = C[[j for j = 1:length(C) if j != pos]]
                (length(lazy_set) > 0) && (lazy_set_g[i] = lazy_set)
            end
        end
    end
    return con_set_g, lazy_set_g
end

function clq_extend_parallel(set::Vector{Vector{Int64}}, cg::SparseMatrixCSC, numthreads::Int64)
    #initialized con_set_g, which is a vector to store extended cliques from threads (save reducing time)
    con_set_g = [Int64[] for _ in 1:length(set)]
    lazy_set_g = Vector{Int64}[]
    m = length(set)
    #width is the lower bound of the number of constraints for each thread
    width = floor(Int64, m / numthreads)
    spl = SpinLock()
    #parallel computation
    @threads for id in 1:numthreads
        let
            #determine indeces of constraints for each thread
            local st = (id - 1) * width + 1
            if id == numthreads
                local en = m
            else
                local en = id * width
            end

            #C is cliques set, and pos is the index, where C[pos] is the longest extended cliques
            C, pos = _extend_unique_clique(set[i], cg)
            #send results back the main core
            con_set_g[i] = C[pos]
            #remove the pos term
            local lazy_set = C[[j for j = 1:length(C) if j != pos]]
            lock(spl)
            (length(lazy_set) > 0) && (append!!(lazy_set_g, lazy_set))
            unlock(spl)
        end
    end
    return con_set_g, lazy_set_g
end

function cliques_extend_serial(set, cg::SparseMatrixCSC)
    return _extend_one_set_I(set, cg)
end


function _extend_one_set_I(set::Vector{Vector{Int64}}, cg::SparseMatrixCSC)
    m = length(set)
    #initialize two sets for cliques added to the original model and the lazy cuts (user cuts)
    con_set = Vector{Int64}[]
    lazy_set = Vector{Int64}[]
    for i in 1:m

        #C is cliques set, and pos is the index, where C[pos] is the longest extended cliques
        C, pos = _extend_unique_clique(set[i], cg)
        push!(con_set, C[pos])
        #remove the pos term
        lazy_set = C[[j for j = 1:length(C) if j != pos]]
    end
    return con_set, lazy_set
end

function _extend_unique_clique(c::AbstractVector, cg::SparseMatrixCSC)
    #degree claculate
    degree_list = [sum(cg[index, :]) + sum(cg[:, index]) for index in c]
    #get the smallest degree variables
    d = c[findmin(degree_list)[2]]
    #initialize an empty list for candidates
    candidate_list = Int64[]

    #find all potential candidates by for loop
    #for all conflict variables to d, there are two part as our represent way of CG
    #first part
    for index in findnz(cg[d, :])[1]
        if !(index in c)
            local candidate_id = true
            #check whether variable index has been connected to all variables from c in CG
            for i in c
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
        if !(index in c)
            local candidate_id = true
            for i in c
                if cg[index, i] + cg[i, index] < 1
                    candidate_id = false
                    break
                end
            end
            (candidate_id) && (push!(candidate_list, index))
        end
    end
    #if there is no candidate, just return the original clique
    (length(candidate_list) > 0) ? (extend_list = [[candidate_list[1]]]) : (return [c], 1)
    #now try to combine candidates to get stronger clique
    for i in candidate_list[2:end]
        #to check wether variable i cannot be combined to any other candidate
        #if yes, we need to generate a clique as original clique + variable i
        local isunique = true
        for (index, k) in enumerate(extend_list)
            #check whether i can be connected to the current clique
            local inthislist = true
            for j in k
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
    end

    #find largest length term in extend_list
    pos = 0
    len = 0
    for (index, k) in enumerate(extend_list)
        (length(k) > len) && ((pos, len) = (index, length(k)))
    end

    return [append!!(k, c) for k in extend_list], pos
end