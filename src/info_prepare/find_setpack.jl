#find_setpack.jl

#find_set_pack() detect cliques and remove them from the original model (cliques not in I)
function find_set_pack(
        con_set::Vector{Vector{Int64}}, con_coef::AbstractVector, 
        con_ub::AbstractVector, con_lb::AbstractVector,
        var_type::AbstractVector, org_to_bin::Dict{Int64, Int64}
    )
    m = length(con_ub)
    #I includes all non-setpack constraints' indices
    I = Int64[]
    """
        J is for the situation like x1 + x2 = 1; 
        then x1 + x2 <= 1 is a clique, but x1 + x2 >= 1 not;
        So we should add x1 + x2 = 1 back to the original model, instead of removing it.
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
                #con_ub[j] = Inf
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
                    #con_lb[j] = -Inf
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
            #for complementary variables
            (coef[i] < 0) && (b += 1)
        end
    end
    (b == 1) ? (return true) : (return false)
end

#represent the clique with variables in binary variables list
function _wrap_set_pack(con::Vector{Int64}, coef::AbstractVector, org_to_bin::Dict{Int64, Int64})
    s = Int64[]
    binary_length = length(org_to_bin)
    for (i, j) in enumerate(con)
        (coef[i] > 0) ? (push!(s, org_to_bin[j])) : (push!(s, org_to_bin[j] + binary_length))
    end
    return s
end
