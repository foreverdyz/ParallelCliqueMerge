#generate_sets.jl

include("detect_cliques.jl")
include("cg_construct.jl")
include("clique_extension.jl")
include("clique_domination.jl")

using SparseArrays

"""
generate_sets() produces cliques sets based on info read from file.
clq_detect_time_limit: time limit for clq extension
cg_length_limit: size limit of clqs in cg construction
cg_term_limit: number of nonzero terms limit in one thread!!! in cg construction
clq_ext_length_limit: size limit of clqs in clq extension
clq_ext_term_limit: number of nonzero terms limit in one thread!!! in clq extension
domain_limit: number of constraints limit in clq domain
"""
function generate_sets(
        set_pack::Vector{Vector{Int64}}, set_pack_new::Vector{Vector{Int64}}, 
        knapsack_set::Vector{Tuple{Any, Any, Any}},
        org_to_bin::Dict{Int64, Int64}, numthreads::Int64,
        clq_detect_time_limit::Real, cg_length_limit::Int64, cg_term_limit::Int64,
        clq_ext_length_limit::Int64, clq_ext_term_limit::Int64, domain_limit::Int64
    )
    #detect clqs from conflicted knapsacks
    if numthreads == 1
        main_cliques, clique_set = detect_cliques(knapsack_set, clq_detect_time_limit); 
    else
        main_cliques, clique_set = detect_cliques_parallel(knapsack_set, numthreads, clq_detect_time_limit);
    end
    #construct conflict graph (cg)
    cg = cg_construct(set_pack, set_pack_new, main_cliques, clique_set, org_to_bin, numthreads, cg_length_limit, cg_term_limit);
    #extend clqs based on cg
    c1, c12 = clique_extension(set_pack, cg, numthreads, clq_ext_length_limit, clq_ext_term_limit)
    c2, c22 = clique_extension(main_cliques, cg, numthreads, clq_ext_length_limit, clq_ext_term_limit)
    c3, c32 = clique_extension(set_pack_new, cg, numthreads, clq_ext_length_limit, clq_ext_term_limit)
    """
        Here we do not extend clqs from clique_set, because the size is ususally too large and similar to c2 and c22. But you can do that with following commented code
        c4, c42 = clique_extension(clique_set, cg, numthreads, clq_ext_length_limit, clq_ext_term_limit)
    """

    """
        Now combing cliques from different sets (from parallel computing)
    """
    iter = length(c12)
    c_temp = Vector{Int64}[]
    for i in 1:iter
        (length(c12[i]) > 0) && (append!!(c_temp, c12[i]))
    end
    c12 = copy(c_temp)
    iter = length(c22)
    c_temp = Vector{Int64}[]
    for i in 1:iter
        (length(c22[i]) > 0) && (append!!(c_temp, c22[i]))
    end
    c22 = copy(c_temp)
    iter = length(c32)
    c_temp = Vector{Int64}[]
    for i in 1:iter
        (length(c32[i]) > 0) && (append!!(c_temp, c32[i]))
    end
    c32 = copy(c_temp)
    """
        iter = length(c42)
        c_temp = Vector{Int64}[]
        for i in 1:iter
            (length(c42[i]) > 0) && (append!!(c_temp, c42[i]))
        end
        c42 = copy(c_temp)
    """
    #clqs domain
    if length(c1) > 1
        if length(c1) < domain_limit
            c1 = clique_domination(c1, numthreads)
        end
    end
    if length(c2) > 1
        if length(c2) < domain_limit
            c2 = clique_domination(c2, numthreads)
        end
    end
    if length(c3) > 1
        if length(c3) < domain_limit
            c3 = clique_domination(c3, numthreads)
        end
    end

    return c1, c2, c3, c4, c12, c22, c32, c42
end
