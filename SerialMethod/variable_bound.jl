#variable_bound.jl

using SparseArrays

function singleton_update(A::AbstractSparseArray, con_ub::AbstractVector, con_lb::AbstractVector, var_ub::AbstractVector, var_lb::AbstractVector, var_type::AbstractVector)
    #get m constraints and n variables
    m, n = size(A)
    #index set for non-singleton constraint
    I = Int64[]
    for j in 1:m
        if length(findnz(A[j, :])[1]) == 1
            i = findnz(A[j, :])[1][1]
            c = A[j, i]
            if c > 0
                (con_ub[j] != Inf) && (var_ub[i] = min(var_ub[i], con_ub[j]/c))
                (con_lb[j] != -Inf) && (var_lb[i] = max(var_lb[i], con_lb[j]/c))
            else
                (con_ub[j] != Inf) && (var_lb[i] = max(var_lb[i], con_ub[j]/c))
                (con_lb[j] != -Inf) && (var_ub[i] = min(var_ub[i], con_lb[j]/c))
            end
            if var_type[i] > 0
                var_ub[i] = floor(Int64, var_ub[i])
                var_lb[i] = ceil(Int64, var_lb[i])
                if var_ub[i] == 1 && var_lb[i] == 0
                    var_type[i] = 2
                end
            end
        else
            push!(I, j)
        end
    end
    return A[I, :], con_ub[I], con_lb[I], var_ub, var_lb, var_type
end

function bound_strength(A::AbstractSparseArray, con_ub::AbstractVector, con_lb::AbstractVector, var_ub::AbstractVector, var_lb::AbstractVector, var_type::AbstractVector)
    #get m constraints and n variables
    m, n = size(A)
    for j in 1:m
        (con_ub[j] != Inf) && (var_ub, var_lb = _unique_ub_con_bound_strength!(A[j, :], con_ub[j], var_ub, var_lb, var_type))
        (con_lb[j] != -Inf) && (var_ub, var_lb = _unique_lb_con_bound_strength!(A[j, :], con_lb[j], var_ub, var_lb, var_type))
    end
    return var_ub, var_lb
end

function _unique_ub_con_bound_strength(a::AbstractVector, b::Real, var_ub::AbstractVector, var_lb::AbstractVector, var_type::AbstractVector)
    delta_list = spzeros(length(var_type))
    for i in findnz(a)[1]
        (a[i] > 0) ? (delta_list[i] = a[i]*var_lb[i]) : (delta_list[i] = a[i]*var_ub[i])
    end
    delta = sum(delta_list)
    for i in findnz(a)[1]
        local x = (b - delta + delta_list[i])/a[i]
        if a[i] > 0
            if var_ub[i] > x
                (var_type[i] > 0) ? (var_ub[i] = floor(Int64, x)) : (var_ub[i] = x)
            end
        else
            if var_lb[i] < x
                (var_type[i] > 0) ? (var_lb[i] = ceil(Int64, x)) : (var_lb[i] = x)
            end
        end
    end
    return var_ub, var_lb
end

function _unique_lb_con_bound_strength(a::AbstractVector, b::Real, var_ub::AbstractVector, var_lb::AbstractVector, var_type::AbstractVector)
    delta_list = spzeros(length(var_type))
    for i in findnz(a)[1]
        (a[i] > 0) ? (delta_list[i] = a[i]*var_ub[i]) : (delta_list[i] = a[i]*var_lb[i])
    end
    delta = sum(delta_list)
    for i in findnz(a)[1]
        local x = (b - delta + delta_list[i])/a[i]
        if a[i] < 0
            if var_ub[i] > x
                (var_type[i] > 0) ? (var_ub[i] = floor(Int64, x)) : (var_ub[i] = x)
            end
        else
            if var_lb[i] < x
                (var_type[i] > 0) ? (var_lb[i] = ceil(Int64, x)) : (var_lb[i] = x)
            end
        end
        return var_ub, var_lb
    end
end