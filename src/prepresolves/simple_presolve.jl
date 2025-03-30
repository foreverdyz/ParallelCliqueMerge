#simple_presolve.jl

using SparseArrays

"""
    simple_presolve(A, con_set, con_coef, con_ub, con_lb, var_ub, var_lb, var_type)

Performs basic presolve techniques on the input problem:
1. Singleton removal: removes constraints with a single variable and strengthens bounds.
2. Bounds strengthening: propagates tightened bounds based on the constraints.
3. Binary re-identification: identifies binary variables and maps them.

# Arguments
- `A::AbstractSparseArray`: The constraint matrix.
- `con_set::Vector{Vector{Int64}}`: List of constraints.
- `con_coef::AbstractVector`: Coefficients of constraints.
- `con_ub::AbstractVector`: Upper bounds of constraints.
- `con_lb::AbstractVector`: Lower bounds of constraints.
- `var_ub::AbstractVector`: Upper bounds of variables.
- `var_lb::AbstractVector`: Lower bounds of variables.
- `var_type::AbstractVector`: Variable types (continuous, integer, binary).

# Returns
Modified constraint matrix, constraint set, constraint coefficients, bounds, and mappings between original and binary variables.
"""
function simple_presolve(A::AbstractSparseArray, con_set::Vector{Vector{Int64}}, con_coef::AbstractVector, con_ub::AbstractVector, con_lb::AbstractVector, var_ub::AbstractVector, var_lb::AbstractVector, var_type::AbstractVector)
    # Step 1: Singleton removal and bounds updating
    I, var_ub, var_lb = _singleton_remove(con_set, con_coef, con_ub, con_lb, var_ub, var_lb)
    
    # Step 2: Strengthen bounds based on the remaining constraints
    con_ub, con_lb, var_ub, var_lb = _bounds_strengthen(con_set, con_coef, con_ub, con_lb, var_ub, var_lb, I, var_type)
    
    # Remove processed constraints
    A = A[I, :]
    con_set = con_set[I]
    con_coef = con_coef[I]
    
    # Step 3: Re-identify binary variables
    var_type, org_to_bin, bin_to_org = _iden_bin_list(var_ub, var_lb, var_type)
    
    return A, con_set, con_coef, con_ub, con_lb, var_ub, var_lb, var_type, org_to_bin, bin_to_org
end

"""
    _singleton_remove(con_set, con_coef, con_ub, con_lb, var_ub, var_lb)

Removes singleton constraints and updates variable bounds accordingly.

# Arguments
- `con_set::Vector{Vector{Int64}}`: List of constraints.
- `con_coef::AbstractVector`: Coefficients of constraints.
- `con_ub::AbstractVector`: Upper bounds of constraints.
- `con_lb::AbstractVector`: Lower bounds of constraints.
- `var_ub::AbstractVector`: Upper bounds of variables.
- `var_lb::AbstractVector`: Lower bounds of variables.

# Returns
- Index set of non-singleton constraints and updated variable bounds.
"""
function _singleton_remove(con_set::Vector{Vector{Int64}}, con_coef::AbstractVector, con_ub::AbstractVector, con_lb::AbstractVector, var_ub::AbstractVector, var_lb::AbstractVector)
    m = length(con_set)
    I = Int64[]  # Index set for non-singleton constraints

    for j in 1:m
        con = con_set[j]
        if length(con) == 1  # Singleton constraint
            i = con[1]
            c = con_coef[j][1]
            
            # Strengthen variable bounds
            if c > 0
                (con_ub[j] != Inf) && (var_ub[i] = min(var_ub[i], con_ub[j]/c))
                (con_lb[j] != -Inf) && (var_lb[i] = max(var_lb[i], con_lb[j]/c))
            else
                (con_ub[j] != Inf) && (var_lb[i] = max(var_lb[i], con_ub[j]/c))
                (con_lb[j] != -Inf) && (var_ub[i] = min(var_ub[i], con_lb[j]/c))
            end
        elseif length(con) > 1
            push!(I, j)
        end
    end
    
    return I, var_ub, var_lb
end

"""
    _bounds_strengthen()
    this function propagates a simple variable bounds strengthening
"""
function _bounds_strengthen(con_set::Vector{Vector{Int64}}, con_coef::AbstractVector, con_ub::AbstractVector, con_lb::AbstractVector, var_ub::AbstractVector, var_lb::AbstractVector, I::Vector{Int64}, var_type::AbstractVector)
    for j in I
        """
            the logic is as follows:
            1, <= and >= will be treated as one inequality
            2, == will be treated as two inequalities separately
        """
        if length(con_set) < 100_000 #avoid exceeding runtime
            (con_ub[j] != Inf) && (var_ub, var_lb = _unique_ub_con_bound_strengthen(con_set[j], con_coef[j], con_ub[j], var_ub, var_lb, var_type))
            (con_lb[j] != -Inf) && (var_ub, var_lb = _unique_ub_con_bound_strengthen(con_set[j], -con_coef[j], -con_lb[j], var_ub, var_lb, var_type))
        end
    end
    #output constraints matrix and bounds with only constraints from I and updated variables bounds
    return con_ub[I], con_lb[I], var_ub, var_lb
end

"""
    _unique_ub_con_bound_strengthen(con, coef, b, var_ub, var_lb, var_type)

Strengthens upper bounds for variables in the constraint `con`.
"""
function _unique_ub_con_bound_strengthen(con::Vector{Int64}, coef::AbstractVector, b::Real, var_ub::AbstractVector, var_lb::AbstractVector, var_type::AbstractVector)
    #ub means we want lower bouond for every variable; we use delta_list to record each lower bound in this constraint
    delta_list = [coef[i] > 0 ? coef[i] * var_lb[con[i]] : coef[i] * var_ub[con[i]] for i in eachindex(con)]
    delta = sum(delta_list)
    (delta == -Inf) && (return var_ub, var_lb)

    for (i, j) in enumerate(con)
        #the new potential bound
        x = (b - delta + delta_list[i]) / coef[i]
        #====
            Here we do not need to udpate the delta, since the updated variable bound is not to do with the delta!
        ====#
        if coef[i] > 0 && var_ub[j] > x + 1e-5
            var_ub[j] = var_type[j] > 0 ? floor(Int64, x) : x
        elseif coef[i] < 0 && var_lb[j] < x - 1e-5
            var_lb[j] = var_type[j] > 0 ? ceil(Int64, x) : x
        end
    end
    return var_ub, var_lb
end

"""
    _iden_bin_list()
    this function implements two operations: first, re-indentify (potential) binary variables from integer variables; second,
build a binary map (including bin to org and org to bin) to help us find a binary variable from original variables and find
the original variable for a binary variable.
"""
function _iden_bin_list(var_ub::AbstractVector, var_lb::AbstractVector, var_type::AbstractVector)
    org_to_bin = Dict{Int64, Int64}()
    bin_to_org = Dict{Int64, Int64}()
    index = 0
     for i in eachindex(var_type)
        if var_type[i] == 1 && var_ub[i] == 1 && var_lb[i] == 0  # Reclassify integer vars as binary
            var_type[i] = 2
            index += 1
            org_to_bin[i] = index
            bin_to_org[index] = i
        elseif var_type[i] == 2
            index += 1
            org_to_bin[i] = index
            bin_to_org[index] = i
        end
    end
    return var_type, org_to_bin, bin_to_org
end