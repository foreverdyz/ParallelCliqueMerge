#process_all_cons.jl

using SparseArrays

#transfer all constraints to pure bianry (with complementary variables) constraints,
#and check whether they are setpack or conflicted knapsack constraints.
function process_all_cons_new(
        con_set::Vector{Vector{Int64}}, con_coef::AbstractVector,
        con_ub::AbstractVector, con_lb::AbstractVector,
        var_ub::AbstractVector, var_lb::AbstractVector,
        var_type::AbstractVector, org_to_bin::Dict{Int64, Int64}
    )
    local set_pack_new = Vector{Int64}[]
    local knapsack_set = Tuple{Any, Any, Any}[]
    m = length(con_set)
    for i in 1:m
        local con, coef = con_set[i], con_coef[i]
        if con_ub[i] != Inf
            #convert to pure bianry and <= constraints
            local bin_con, bin_coef, bin_rhs = _convert_to_standard(
                con, coef, con_ub[i],
                var_ub, var_lb,
                var_type, org_to_bin
            )
            if length(bin_con) > 1
                #check whether it is a clique or a conflicted knapsack (in which we can detect new clqs later)
                local id = _check_conflict(bin_coef, bin_rhs)
                if id == 2
                    push!(set_pack_new, bin_con)
                elseif id == 1
                    push!(knapsack_set, (bin_con, bin_coef, bin_rhs))
                end
            end
        end
        
        if con_lb[i] != -Inf
            local bin_con, bin_coef, bin_rhs = _convert_to_standard(
                con, -coef, -con_lb[i],
                var_ub, var_lb,
                var_type, org_to_bin
            )
            if length(bin_con) > 1
                local id = _check_conflict(bin_coef, bin_rhs)
                if id == 2
                    push!(set_pack_new, bin_con)
                elseif id == 1
                    push!(knapsack_set, (bin_con, bin_coef, bin_rhs))
                end
            end
        end
    end
    return set_pack_new, knapsack_set
end

#next step: convert to standard; a constraint, remove all non binary term, and replace negative coefficients
#term to their complementary variables
#this function is only for ax <= b, so please transfer ax >= b to -ax <= -b
function _convert_to_standard(
        con::Vector{Int64}, coef::AbstractVector, b::Real,
        var_ub::AbstractVector, var_lb::AbstractVector,
        var_type::AbstractVector, org_to_bin::Dict{Int64, Int64}
    )
    binary_length = length(org_to_bin)
    #initilize a sparse vector for the pure binary constraint (with complementary binary variables)
    #bin_con = spzeros(2*binary_length)
    bin_con = Int64[]
    bin_coef = []
    #initiliaze a new rhs of the constraint
    bin_rhs = copy(b)
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

#output 0: nothing happens, 1: conflicted knapsack constraint, 2: set pack
function _check_conflict(coef::AbstractVector, b::Real)
    s = 0
    for i in coef
        s += min(0, i)
    end
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
        #if i > largest_two_coef[2]
        #    (i > largest_two_coef[1]) ? (largest_two_coef = [i, largest_two_coef[1]]) : (largest_two_coef[2] = i)
        #end
    end
    (is_setpack) && (return 2)
    (sum(partialsort(coef, rev = true, 1:2)) > b + 0.0001) ? (return 1) : (return 0)
end