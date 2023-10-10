#find_set_pack.jl

using SparseArrays

function find_set_pack(A::AbstractSparseArray, con_ub::AbstractVector, con_lb::AbstractVector, var_type::AbstractVector)
    m = size(A)[1]
    I = Int64[]
    S = Int64[]
    for j in 1:m
        is_set_pack = true
        (con_ub[j] != 1) && (is_set_pack = false)
        if is_set_pack
            (con_lb[j] != -Inf) && (is_set_pack = false)
        end
        if is_set_pack
            for i in findnz(A[j,:])[1]
                if var_type[i] != 2
                    is_set_pack = false
                    break
                end
                if abs(A[j, i]) != 1
                    is_set_pack = false
                    break
                end
            end 
        end
        (is_set_pack) ? (push!(S, j)) : (push!(I, j))
    end
    return S
end