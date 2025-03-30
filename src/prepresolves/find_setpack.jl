#find_setpack.jl

"""
    find_setpack()
Find original setpacking constraints that can be replaced by extended cliques later.
"""
function find_setpack(
        con_set::Vector{Vector{Int64}}, con_coef::AbstractVector, 
        con_ub::AbstractVector, con_lb::AbstractVector,
        var_type::AbstractVector, org_to_bin::Dict{Int64, Int64}
    )
    m = length(con_ub)
    #I includes all non-setpack constraints' indices
    I, J = Int64[], Int64[]
    #set_pack includes bin numbers of terms in setpack constraints
    set_pack = Vector{Int64}[]
    for j in 1:m
        has_setpack = false
        if con_ub[j] != Inf
            if _check_set_pack(con_set[j], con_coef[j], con_ub[j], var_type)
                push!(set_pack, _wrap_set_pack(con_set[j], con_coef[j], org_to_bin))
                has_setpack = true
            else
                push!(I, j)
            end
        end
        if con_lb[j] != -Inf
             if has_setpack
                push!(J, j) #here is to avoid some setpacking constraints considered again in later step
            elseif _check_set_pack(con_set[j], -con_coef[j], -con_lb[j], var_type)
                push!(set_pack, _wrap_set_pack(con_set[j], -con_coef[j], org_to_bin))
            else
                # Avoid pushing j into I twice
                if isempty(I) || I[end] != j
                    push!(I, j)
                end
            end
        end
    end
    
    return I, J, set_pack, con_ub, con_lb
end

"""
    _check_set_pack()
Checks if a constraint qualifies as a set-packing constraint.
"""
function _check_set_pack(con::Vector{Int64}, coef::AbstractVector, b::Real, var_type::AbstractVector)
    for (i, j) in enumerate(con)
        if var_type[j] < 2 || abs(coef[i]) != 1
            return false
        end
        (coef[i] < 0) && (b += 1)
    end
    return b == 1
end

"""
    _wrap_set_pack()
Transforms the set-packing constraint into a set of binary variables.
"""
function _wrap_set_pack(con::Vector{Int64}, coef::AbstractVector, org_to_bin::Dict{Int64, Int64})
    s = Int64[]
    binary_length = length(org_to_bin)
    for (i, j) in enumerate(con)
        (coef[i] > 0) ? (push!(s, org_to_bin[j])) : (push!(s, org_to_bin[j] + binary_length))
    end
    return s
end