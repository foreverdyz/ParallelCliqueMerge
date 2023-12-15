#remap_cons.jl

#before using cliques for the original MIP, we need to remap these binary variables back to their orignal positions
function remap_cons(set::Vector{Vector{Int64}}, bin_to_org::Dict{Int64, Int64}, var_lb::AbstractVector)
    m = length(set)
    #initialize a set to store remapped constraints
    gen_set = []
    binary_number = length(bin_to_org)
    for j in 1:m
        c = spzeros(length(var_lb))
        for i in set[j]
            #i >  binary_number => i is a complementary variable
            (i <= binary_number) ? (c[bin_to_org[i]] = 1) : (c[bin_to_org[i - binary_number]]  = -1)
        end
        push!(gen_set, c)
    end
    return gen_set
end