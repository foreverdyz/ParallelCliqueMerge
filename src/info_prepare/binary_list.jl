#binary_list.jl

"""
    binary_list() build two dictionaries:
    one for original varibles to binary variables,
    and one for binary variables to original variables
"""
function binary_list(var_type::AbstractVector)
    org_to_bin = Dict{Int64, Int64}()
    bin_to_org = Dict{Int64, Int64}()
    index = 0
    for (i, v) in enumerate(var_type)
        if v == 2
            index += 1
            org_to_bin[i] = index
            bin_to_org[index] = i
        end
    end
    return org_to_bin, bin_to_org
end