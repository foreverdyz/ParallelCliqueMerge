#simple_presolve.jl

using SparseArrays

"""
The simple presolve includes one round of singleton removing, one-row variables bounds strengthenin, and binary re-identification
"""
function simple_presolve(
    A::AbstractSparseArray, con_set::Vector{Vector{Int64}}, con_coef::AbstractVector, 
    con_ub::AbstractVector, con_lb::AbstractVector, var_ub::AbstractVector, var_lb::AbstractVector, 
    var_type::AbstractVector, max_tolerance::Float64
    )
    I, var_ub, var_lb, var_type = _singleton_remove(con_set, con_coef, con_ub, con_lb, var_ub, var_lb, var_type)

    con_ub, con_lb, var_ub, var_lb = _bounds_strengthen(
        con_set, con_coef, con_ub, con_lb, 
        var_ub, var_lb, I, 
        var_type, max_tolerance
        )
    
    #I incldues indeces of all remaining constraints
    A = A[I, :]
    con_set = con_set[I]
    con_coef = con_coef[I]

    #indentify more binary variables from integer variables
    for i in 1:length(var_ub)
        if var_type[i] == 1
            if var_ub[i] == 1 && var_lb[i] == 0
                var_type[i] = 2
            end
        end
    end
    return A, con_set, con_coef, con_ub, con_lb, var_ub, var_lb, var_type
end

"""
    singleton_remove()
    this function removes singleton from the constraints matrix and update variable bounds 
"""
function _singleton_remove(con_set::Vector{Vector{Int64}}, con_coef::AbstractVector, con_ub::AbstractVector, con_lb::AbstractVector, var_ub::AbstractVector, var_lb::AbstractVector, var_type::AbstractVector)
    #get m constraints
    m = length(con_set)
    #index set for non-singleton constraint
    I = Int64[]
    #check all constraints
    for j in 1:m
        local con = con_set[j]
        if length(con) == 1 #which implies len(c) is 1
            #i is the index of the only variable in the j-th constraints
            local i = con[1]
            #get the coefficient
            local c = con_coef[j][1]
            #then we can strengthen variables' bounds
            if c > 0
                (con_ub[j] != Inf) && (var_ub[i] = min(var_ub[i], con_ub[j]/c))
                (con_lb[j] != -Inf) && (var_lb[i] = max(var_lb[i], con_lb[j]/c))
            else
                (con_ub[j] != Inf) && (var_lb[i] = max(var_lb[i], con_ub[j]/c))
                (con_lb[j] != -Inf) && (var_ub[i] = min(var_ub[i], con_lb[j]/c))
            end
            #if this variable is integer (1 or 2), we can round the upper and lower bounds
            if var_type[i] > 0
                (var_ub[i] != Inf) && (var_ub[i] = floor(Int64, var_ub[i]))
                (var_lb[i] != -Inf) && (var_lb[i] = ceil(Int64, var_lb[i]))
                #some integer variables can be noted as binary variable
                if var_ub[i] == 1 && var_lb[i] == 0
                    var_type[i] = 2
                end
            end
        elseif length(con) > 1 #remove empty constraints
            #j is not a singleton constraint; thus push it to the index set I
            push!(I, j)
        end
    end
    #output the index set I and updated varaibles' bounds
    return I, var_ub, var_lb, var_type
end

"""
    bounds_stren()
    this function propagates a simple variable bounds strengthening
"""
function _bounds_strengthen(
    con_set::Vector{Vector{Int64}}, con_coef::AbstractVector, con_ub::AbstractVector, con_lb::AbstractVector, 
    var_ub::AbstractVector, var_lb::AbstractVector, I::Vector{Int64}, 
    var_type::AbstractVector, max_tolerance::Float64
    )
    for j in I
        """
            the logic is as follows:
            1, <= and >= will be treated as one inequality
            2, == will be treated as two inequalities separately
        """
        (con_ub[j] != Inf) && (var_ub, var_lb = _unique_ub_con_bound_strengthen(con_set[j], con_coef[j], con_ub[j], var_ub, var_lb, var_type, max_tolerance))
        (con_lb[j] != -Inf) && (var_ub, var_lb = _unique_ub_con_bound_strengthen(con_set[j], -con_coef[j], -con_lb[j], var_ub, var_lb, var_type, max_tolerance))
    end
    #output constraints matrix and bounds with only constraints from I and updated variables bounds
    return con_ub[I], con_lb[I], var_ub, var_lb
end

function _unique_ub_con_bound_strengthen(
    con::Vector{Int64}, coef::AbstractVector, 
    b::Real,var_ub::AbstractVector, var_lb::AbstractVector, 
    var_type::AbstractVector, max_tolerance::Float64
    )
    #ub means we want lower bouond for every variable; we use delta_list to record each lower bound in this constraint
    delta_list = zeros(length(con))
    for (i, j) in enumerate(con)
        #calculate lower bound for each variables
        (coef[i] > 0) ? (delta_list[i] = coef[i]*var_lb[j]) : (delta_list[i] = coef[i]*var_ub[j])
    end
    delta = sum(delta_list)
    (delta == -Inf) && (return var_ub, var_lb)
    for (i, j) in enumerate(con)
        #the new potential bound
        local x = (b - delta + delta_list[i])/coef[i]
        #====
            Here we do not need to udpate the delta, since the updated variable bound is not to do with the delta!
        ====#
        if coef[i] > 0
            if var_ub[j] > x + max_tolerance
                (var_type[j] > 0) ? (var_ub[j] = floor(Int64, x)) : (var_ub[j] = x)
            end
        else
            if var_lb[j] < x - max_tolerance
                (var_type[j] > 0) ? (var_lb[j] = ceil(Int64, x)) : (var_lb[j] = x)
            end
        end
    end
    #output
    return var_ub, var_lb
end