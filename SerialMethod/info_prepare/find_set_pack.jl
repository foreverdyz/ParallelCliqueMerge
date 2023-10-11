#find_set_pack.jl

using SparseArrays

function find_set_pack(A::AbstractSparseArray, con_ub::AbstractVector, con_lb::AbstractVector, var_type::AbstractVector, org_to_bin::Dict{Int64,Int64})
    m = size(A)[1]
    I = Int64[]
    set_pack = Vector{Int64}[]
    for j in 1:m
        if _check_set_pack(A[j, :], con_ub[j], con_lb[j], var_type)
            push!(set_pack, _wrap_set_pack(A[j, :], org_to_bin))
        else
            push!(I, j)
        end
    end
    return A[I, :], con_ub[I], con_lb[I], set_pack
end

function _check_set_pack(a::AbstractVector, ub::Real, lb::Real, var_type::AbstractVector)
    (ub != 1) && (return false)
    (lb != -Inf) && (return false)
    for i in findnz(a)[1]
        if var_type[i] != 2
            return false
        end
        if abs(a[i]) != 1
            return false
        end
    end
    return true
end

function _wrap_set_pack(a::AbstractVector, org_to_bin::Dict{Int64,Int64})
    set_pack = Int64[]
    for i in findnz(a)[1]
        if a[i] == 1
            push!(set_pack, org_to_bin[i])
        else
            push!(set_pack, 2*org_to_bin[i])
        end
    end
    return set_pack
end