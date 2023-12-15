#find_set_pack.jl

using SparseArrays

#this funciton tries to find all set-pack constraints
function find_set_pack(A::AbstractSparseArray, con_ub::AbstractVector, con_lb::AbstractVector, var_type::AbstractVector, org_to_bin::Dict{Int64,Int64})
    #get the number of constraints
    m = size(A)[1]
    #initialize an index set to cache non set-pack constraints
    I = Int64[]
    #initialize a set to cache set-pack constraints
    set_pack = Vector{Int64}[]
    #for loop to check all constraints
    for j in 1:m
        #check whether constraint j is set-pack or not
        if _check_set_pack(A[j, :], con_ub[j], con_lb[j], var_type)
            #if j is a set-pack, transfer it to a pure binary form and store to set_pack set
            push!(set_pack, _wrap_set_pack(A[j, :], org_to_bin))
        else
            #if j is not a set-pack, cache it to the index set I
            push!(I, j)
        end
    end
    #output non set-pack constraints and bonds, and set_pack
    return A[I, :], con_ub[I], con_lb[I], set_pack
end

#this functions tries to identify set-pack constraint
function _check_set_pack(a::AbstractVector, ub::Real, lb::Real, var_type::AbstractVector)
    #set-pack has form like <= 1
    (ub != 1) && (return false)
    #avoid ax = 1, since it is stronger than ax <= 1; thus we cannot remove it from the original problem
    (lb != -Inf) && (return false)
    #check all coefficients
    for i in findnz(a)[1]
        #if there are non-binary variables, false
        if var_type[i] != 2
            return false
        end
        #if the coefficient is not 1, false
        if a[i] != 1
            return false
        end
    end
    return true
end

#this function tries to transfer a set-pack to
function _wrap_set_pack(a::AbstractVector, org_to_bin::Dict{Int64,Int64})
    #initialize an empty set to store indeces of the set-pack
    set_pack = Int64[]
    for i in findnz(a)[1]
        #for i, map it to pure binary list
        push!(set_pack, org_to_bin[i])
    end
    return set_pack
end