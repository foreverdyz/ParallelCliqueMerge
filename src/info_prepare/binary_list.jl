#binary_list.jl

#this function will build two dictionaries for map between the original variable list and a pure binary list
function binary_list(var_type::AbstractVector)
    #initialize two empty dictionaries
    org_to_bin = Dict{Int64,Int64}()
    bin_to_org = Dict{Int64,Int64}()
    #index is for the number of variable in the pure binary list
    index = 0
    for (i, v) in enumerate(var_type)
        if v == 2
            index += 1
            org_to_bin[i] = index
            bin_to_org[index] = i
        end
    end
    #output two dictionaries
    return org_to_bin, bin_to_org
end