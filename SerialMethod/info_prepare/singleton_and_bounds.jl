#singleton_and_bounds.jl

using SparseArrays

function singleton_and_bounds(A::AbstractSparseArray, con_ub::AbstractVector, con_lb::AbstractVector, var_ub::AbstractVector, var_lb::AbstractVector, var_type::AbstractVector)
    #get m constraints
    m, n = size(A)
    #index set for non-singleton constraint
    I = Int64[]
    for j in 1:m
        #singleton implies there is only one varaible in this constraint
        if length(findnz(A[j, :])[1]) == 1
            #i is the only variable
            i = findnz(A[j, :])[1][1]
            #get the coefficient
            c = A[j, i]
            #we can strengthen variables' bounds
            if c > 0
                (con_ub[j] != Inf) && (var_ub[i] = min(var_ub[i], con_ub[j]/c))
                (con_lb[j] != -Inf) && (var_lb[i] = max(var_lb[i], con_lb[j]/c))
            else
                (con_ub[j] != Inf) && (var_lb[i] = max(var_lb[i], con_ub[j]/c))
                (con_lb[j] != -Inf) && (var_ub[i] = min(var_ub[i], con_lb[j]/c))
            end
            #if this variable is integer (1 or 2), we can round the upper and lower bounds
            if var_type[i] > 0
                #some integer variables can be noted as binary variable
                if var_ub[i] == 1 && var_lb[i] == 0
                    var_type[i] = 2
                end
            end
        else
            #j is not a singleton constraint; thus push it to the index set I
            push!(I, j)
        end
    end
    for j in I
        (con_ub[j] != Inf) && (var_ub, var_lb = _unique_ub_con_bound_strengthen(A[j, :], con_ub[j], var_ub, var_lb, var_type))
        (con_lb[j] != -Inf) && (var_ub, var_lb = _unique_lb_con_bound_strengthen(A[j, :], con_lb[j], var_ub, var_lb, var_type))
    end
    #update varaible type again
    for i in 1:n
        if var_type[i] > 0
            var_ub[i] = floor(Int64, var_ub[i])
            var_lb[i] = ceil(Int64, var_lb[i])
            #some integer variables can be noted as binary variable
            if var_ub[i] == 1 && var_lb[i] == 0
                var_type[i] = 2
            end
        end
    end
    #return the new constraint matrix with only I, its corresponding constraints' upper and lower bounds, and new varaibles's bounds and types
    return A[I, :], con_ub[I], con_lb[I], var_ub, var_lb, var_type
end

function _unique_ub_con_bound_strengthen(a::AbstractVector, b::Real, var_ub::AbstractVector, var_lb::AbstractVector, var_type::AbstractVector)
    #ub means we want lower bouond for every variable; we use delta_list to record each lower bound in this constraint
    delta_list = spzeros(length(var_type))
    for i in findnz(a)[1]
        #calculate lower bound for each variables
        (a[i] > 0) ? (delta_list[i] = a[i]*var_lb[i]) : (delta_list[i] = a[i]*var_ub[i])
    end
    delta = sum(delta_list)
    for i in findnz(a)[1]
        #new potential bounds
        local x = (b - delta + delta_list[i])/a[i]
        if a[i] > 0
            #note that, we will never change delta_list in one constraint, very magical
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

function _unique_lb_con_bound_strengthen(a::AbstractVector, b::Real, var_ub::AbstractVector, var_lb::AbstractVector, var_type::AbstractVector)
    delta_list = spzeros(length(var_type))
    for i in findnz(a)[1]
        #calculate upper bound for each variable
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