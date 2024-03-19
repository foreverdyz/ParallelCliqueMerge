#generate_sets.jl

include("detect_cliques.jl")
include("cg_construct.jl")
include("clique_extension.jl")
include("clique_domination.jl")

using SparseArrays

"""
    generate_sets()
    generates all clique sets we want to build the reduced model
"""
function generate_sets(
        set_pack::Vector{Vector{Int64}}, set_pack_new::Vector{Vector{Int64}}, 
        knapsack_set::Vector{Tuple{Any, Any, Any}},
        org_to_bin::Dict{Int64, Int64}, numthreads::Int64,
        max_length::Int64, max_nonzero::Int64, max_domian::Int64
    )
    #detect cliques from conflicting knapsacks
    main_cliques, clique_set = detect_cliques_parallel(knapsack_set, numthreads);
    #construct conflict graph
    cg = cg_construct(set_pack, set_pack_new, main_cliques, clique_set, org_to_bin, numthreads, max_length);
    #extend cliques
    c1, c12 = clique_extension(set_pack, cg, numthreads, max_nonzero)
    c2, c22 = clique_extension(main_cliques, cg, numthreads, max_nonzero)
    c3, c32 = clique_extension(set_pack_new, cg, numthreads, max_nonzero)
    #do not extend c4
    c4 = clique_set;
    c42 = Vector{Int64}[]; 

    if length(c1) > 1
        if length(c1) < max_domian
            c1 = clique_domination(c1, numthreads)
        end
    end
    if length(c2) > 1
        if length(c2) < max_domian
            c2 = clique_domination(c2, numthreads)
        end
    end
    if length(c3) > 1
        if length(c3) < max_domian
            c3 = clique_domination(c3, numthreads)
        end
    end

    return c1, c2, c3, c4, c12, c22, c32, c42
end