#process_all_cons.jl

using SparseArrays

"""
    process_all_cons()
Transforms all constraints into pure binary constraints (with complementary variables) and checks whether they are set-packing or conflicted knapsack constraints.
"""
function process_all_cons(
        con_set::Vector{Vector{Int64}}, con_coef::AbstractVector, con_ub::AbstractVector, con_lb::AbstractVector,
        var_ub::AbstractVector, var_lb::AbstractVector, var_type::AbstractVector, org_to_bin::Dict{Int64, Int64}
    )
    set_pack_new, knapsack_set = Vector{Int64}[], Tuple{Any, Any, Any}[]
    m = length(con_set)
    for i in 1:m
        if length(con_set[i]) < 100_000
            con, coef = con_set[i], con_coef[i]
            if con_ub[i] != Inf
                bin_con, bin_coef, bin_rhs = _convert_to_standard(
                    con, coef, con_ub[i], var_ub, var_lb, var_type, org_to_bin
                )
                if length(bin_con) > 1
                    id = _check_conflict(bin_coef, bin_rhs)
                    if id == 2
                        push!(set_pack_new, bin_con) 
                    elseif id == 1
                        push!(knapsack_set, (bin_con, bin_coef, bin_rhs))
                    end
                end
            end

            if con_lb[i] != -Inf
                bin_con, bin_coef, bin_rhs = _convert_to_standard(
                    con, -coef, -con_lb[i], var_ub, var_lb, var_type, org_to_bin
                )
                if length(bin_con) > 1
                    id = _check_conflict(bin_coef, bin_rhs)
                    if id == 2
                        push!(set_pack_new, bin_con)
                    elseif id == 1
                        push!(knapsack_set, (bin_con, bin_coef, bin_rhs))
                    end
                end
            end
        end
    end

    return set_pack_new, knapsack_set
end

"""
    _convert_to_standard()
Converts a constraint to its standard form by removing non-binary terms and replacing negative coefficients with their complementary variables.
This function is for constraints of the form `Ax ≤ b`.
"""
function _convert_to_standard(
        con::Vector{Int64}, coef::AbstractVector, b::Real,
        var_ub::AbstractVector, var_lb::AbstractVector,
        var_type::AbstractVector, org_to_bin::Dict{Int64, Int64}
    )
    binary_length = length(org_to_bin)
    bin_con, bin_coef = Int64[], Real[]
    bin_rhs = b  # Initialize new RHS
    
    #for all nonzero terms in con
    for (i, j) in enumerate(con)
        #if j is not a binary variable, move it to rhs
        if var_type[j] < 2
            if coef[i] > 0
                (var_lb[j] == -Inf) ? (return Int64[], Int64[], 0) : (bin_rhs -= coef[i] * var_lb[j])
            else
                (var_ub[j] == Inf) ? (return Int64[], Int64[], 0) : (bin_rhs -= coef[i] * var_ub[j])
            end
        else
            if coef[i] > 0 
                push!(bin_con, org_to_bin[j])
                push!(bin_coef, coef[i])
            else
                push!(bin_con, org_to_bin[j] + binary_length)
                push!(bin_coef, -coef[i])
                bin_rhs -= coef[i]
            end
        end
    end
    return bin_con, bin_coef, bin_rhs
end

"""
    _check_conflict()
Classifies a constraint as either set-packing or a conflicted knapsack constraint.
Returns:
- 0: no conflict
- 1: conflicted knapsack constraint
- 2: set-packing constraint
"""
function _check_conflict(coef::AbstractVector, b::Real)
    s = sum(min(0, c) for c in coef)  # Sum of negative coefficients
    #s > b, infeasible; s = b, implies variables should be fixed 
    (s >= b) && (return 0)
    is_setpack = true
    (b != 1) && (is_setpack = false)
    #largest_two_coef includes two elements, term1 and term2, where term1 >  term2
    #largest_two_coef = [-Inf, -Inf]
    for i in coef
        #check whether this constraint is a setpack
        if is_setpack
            (i != 1) && (is_setpack = false)
        end
    end
    (is_setpack) && (return 2)
    (sum(partialsort(coef, rev = true, 1:2)) > b + 1e-5) ? (return 1) : (return 0)
end