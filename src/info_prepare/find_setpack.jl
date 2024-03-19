#find_setpack.jl

"""
    find_set_pack()
    Find org set pack from MIP
"""
function find_set_pack(
        con_set::Vector{Vector{Int64}}, con_coef::AbstractVector, 
        con_ub::AbstractVector, con_lb::AbstractVector,
        var_type::AbstractVector, org_to_bin::Dict{Int64, Int64}
    )
    m = length(con_ub)
    #I includes all non-setpack constraints' indices
    I = Int64[]
    """
        J is for these equations. For example:
            x1 + x2 == 1 can generate a set pack x1 + x2 <=1;
        However, we cannot replace x1 + x2 == 1 by x1 + x2 <= 1;
        and we also do not want to find implied set pack from x1 + x2 == 1 (repeat);
        thus we add x1 + x2 == 1 to J, 
        which will not be processed again but will be added to the final model.
    """
    J = Int64[]
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
                push!(J, j)
            else
                if _check_set_pack(con_set[j], -con_coef[j], -con_lb[j], var_type)
                    push!(set_pack, _wrap_set_pack(con_set[j], -con_coef[j], org_to_bin))
                else
                    #avoid push j to I twice
                    if length(I) > 0
                        (I[end] != j) && (push!(I, j))
                    else
                        push!(I, j)
                    end
                end
            end
        end
    end
    return I, J, set_pack, con_ub, con_lb
end

function _check_set_pack(con::Vector{Int64}, coef::AbstractVector, b::Real, var_type::AbstractVector)
    for (i, j) in enumerate(con)
        if var_type[j] < 2
            return false
        else
            (abs(coef[i]) != 1) && (return false)
            (coef[i] < 0) && (b += 1)
        end
    end
    (b == 1) ? (return true) : (return false)
end

function _wrap_set_pack(con::Vector{Int64}, coef::AbstractVector, org_to_bin::Dict{Int64, Int64})
    s = Int64[]
    binary_length = length(org_to_bin)
    for (i, j) in enumerate(con)
        (coef[i] > 0) ? (push!(s, org_to_bin[j])) : (push!(s, org_to_bin[j] + binary_length))
    end
    return s
end
