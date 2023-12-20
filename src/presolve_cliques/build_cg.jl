#build_cg.jl

using Base.Threads
using SparseArrays

"""
Here we need number of binary variables
"""

#build a conflict graph (cg) with only 1 thread (serial)
function build_cg_serial(set_pack::Vector{Vector{Int64}}, binary_number::Int64)
    #initialize a cg
    cg = _init_cg(binary_number)

    m = length(set_pack)
    #build cg based on all set-packing constraints
    for j in 1:m
        #check all pairs of elements in a constraint
        for i in 2:length(set_pack[j])
            #index-1 to avoid diag terms
            for k in 1:(i-1)
                #cg is a triangle matrix
                (set_pack[j][i] > set_pack[j][k]) ? (cg[set_pack[j][i], set_pack[j][k]] = 1) : (cg[set_pack[j][k], set_pack[j][i]] = 1)
            end
        end
    end
    return cg
end

#update the conflict graph (cg) with 1 thread (serial)
function update_cg_serial(set::Vector{Vector{Int64}}, binary_number::Int64, cg::SparseMatrixCSC)
    m = length(set)
    #update cg based on all set-packing (set) constraints
    for j in 1:m
        for i in 2:length(set[j])
            #index-1 to avoid diag terms
            for k in 1:(i-1)
                (set[j][i] > set[j][k]) ? (cg[set[j][i], set[j][k]] = 1) : (cg[set[j][k], set[j][i]] = 1)
            end
        end
    end
    return cg
end

#build a conflict graph (cg) with multiple threads (distributed structure)
function build_cg_distributed(set_pack::Vector{Vector{Int64}}, binary_number::Int64, numthreads::Int64)
    #initialize a cg
    cg_init = _init_cg(binary_number)

    m = length(set_pack)
    #if m = 0, scheduling a parallel computation will waste some runtime
    (m < 1) && (return cg_init)
    #m is the upper bound of the number of threads
    numthreads = min(m, numthreads)
    #the lower bound of the number of constraint for one thread
    width = floor(Int64, m / numthreads)

    #we pre-build a vector of cgs to store cgs from different threads (save lock time in reducing)
    cg_unit = spzeros(Int64, 2 * binary_number, 2 * binary_number) #this is an empty cg
    cg_g = [cg_unit for _ in 1:numthreads]
    #parallel computation
    @threads for id in 1:numthreads
        let
            #determine the indeces of constraints for this thread
            local st = (id - 1) * width + 1
            if id == numthreads
                local en = m
            else
                local en = id * width
            end

            #since in each thread, we start from an empty sparse matrix,
            #we can build the sparse matrix by this way:
            #a - vector for indeces, b - vector for indeces, and a unit vector for values
            local a = Int64[]
            local b = Int64[]
            #s represents the total number of elements
            local s = 0
            for j in st:en
                #here we want diag element keeps zero
                for i in 2:length(set_pack[j])
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
            end
            #here we add an element like: cg[2*binary_number, 2*binary_number] = 0
            #to guarantee the size of the sparse matrix is (2*binary_number, 2*binary_number)
            push!(a, 2 * binary_number)
            push!(b, 2 * binary_number)
            c = ones(Int64, s)
            push!(c, 0)
            #build cg (sparse(a, b, c)), remove zeors, and reduce to the main core
            cg_g[id] = dropzeros!(sparse(a, b, c))
        end
    end
    #1 thread => does not need to combine cgs by binary combination
    (numthreads == 1) && (return cg_g[1] .| cg_init)
    #start the binary combination
    #L is the number of levels
    L = ceil(Int64, log2(numthreads))
    K = numthreads
    for _ in 1:L
        local k = ceil(Int64, K / 2)
        @threads for i in 1:k
            #combine two cgs
            (i + k <= K) && (cg_g[i] = cg_g[i] .| cg_g[i+k])
        end
        K = k
    end
    return cg_g[1] .| cg_init
end

#update a conflict graph (cg) with multiple threads (distributed structure)
function update_cg_distributed(set::Vector{Vector{Int64}}, binary_number::Int64, cg_init::SparseMatrixCSC, numthreads::Int64)
    m = length(set)
    #avoid scheduling parallel computation for no costraints situation
    (m < 1) && (return cg_init)
    #m is the upper bound of the number of threads
    numthreads = min(m, numthreads)
    #width is the lower bound of the number of constraint for one thread
    width = floor(Int64, m / numthreads)

    #we pre-build a vector of cgs to store cgs from different threads (save lock time in reducing)
    cg_unit = spzeros(Int64, 2 * binary_number, 2 * binary_number) #an empty cg
    cg_g = [cg_unit for _ in 1:numthreads]

    #parallel computation
    @threads for id in 1:numthreads
        let
            #determine indeces for constraints in this thread
            local st = (id - 1) * width + 1
            if id == numthreads
                local en = m
            else
                local en = id * width
            end

            #build cg by a, b (vectors of indeces), and a unit vector
            local a = Int64[]
            local b = Int64[]
            #s represents the total number of elements
            local s = 0
            for j in st:en
                #here we want diag element keeps zero
                for i in 2:length(set[j])
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
            end
            #here we add an element like: cg[2*binary_number, 2*binary_number] = 0
            #to guarantee the size of the sparse matrix is (2*binary_number, 2*binary_number)
            push!(a, 2 * binary_number)
            push!(b, 2 * binary_number)
            c = ones(Int64, s)
            push!(c, 0)
            #build cg (sparse(a, b, c)), remove zeors, and reduce to the main core
            cg_g[id] = dropzeros!(sparse(a, b, c))
        end
    end

    #1 thread => does not need to combine cgs by binary combination
    (numthreads == 1) && (return cg_g[1] .| cg_init)
    #start the binary combination
    #L is the number of levels
    L = ceil(Int64, log2(numthreads))
    K = numthreads
    for _ in 1:L
        local k = ceil(Int64, K / 2)
        @threads for i in 1:k
            #combine two cgs
            (i + k <= K) && (cg_g[i] = cg_g[i] .| cg_g[i+k])
        end
        K = k
    end
    return cg_g[1] .| cg_init
end

#initialize a cg which has included conflicts between the variable and its complementary variable
function _init_cg(binary_number::Int64)
    #initialize an empty sparse matrix
    cg = spzeros(Int64, 2 * binary_number, 2 * binary_number)
    #x + bar(x) = 1, so x + bar(x) <= 1
    for i in 1:binary_number
        cg[i+binary_number, i] = 1
    end
    return cg
end
