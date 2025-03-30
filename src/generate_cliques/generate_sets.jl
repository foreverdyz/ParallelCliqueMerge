#generate_sets.jl

include("detect_cliques.jl");
include("cg_construct.jl");
include("clique_extension.jl");
include("clique_domination.jl");

using Random

function generate_sets(
        set_pack::Vector{Vector{Int64}}, set_pack_new::Vector{Vector{Int64}}, 
        knapsack_set::Vector{Tuple{Any, Any, Any}}, org_to_bin::Dict{Int64, Int64}, numthreads::Int64, 
        knapsack_size::Int64, clique_size::Int64, process_limit_cg::Int64,  time_limit_cg::Real, 
        process_limit_extend::Int64, time_limit_extend::Real
    )
    Random.seed!(12)
    main_cliques, clique_set = detect_cliques(knapsack_set, numthreads, knapsack_size);
    cg, cg_trans = cg_construct(
        set_pack, set_pack_new, main_cliques, clique_set, org_to_bin,
        numthreads, clique_size, process_limit_cg, time_limit_cg
    );
    time_limit = time_limit_extend/log2(numthreads+1);
    t = time();
    c1, c12 = clique_extension(set_pack, cg, cg_trans, numthreads,  clique_size, process_limit_extend, time_limit);
    t_used = time() - t;
    if time_limit - t_used > 0
        c2, c22 = clique_extension(main_cliques, cg, cg_trans, numthreads,  clique_size, process_limit_extend, time_limit - t_used);
    else
        c2, c22 = main_cliques, Vector{Int64}[];
    end
    t_used = time() - t;
    if time_limit - t_used > 0
        c3, c32 = clique_extension(set_pack_new, cg, cg_trans, numthreads,  clique_size, process_limit_extend, time_limit - t_used);
    else
        c3, c32 = set_pack_new, Vector{Int64}[];
    end
    t_used = time() - t;
    if time_limit - t_used > 0
        c4, c42 = clique_extension(clique_set, cg, cg_trans, numthreads,  clique_size, process_limit_extend, time_limit - t_used);
    else
        c4, c42 = clique_set, Vector{Int64}[];
    end
    
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
    iter = length(c42)
    c_temp = Vector{Int64}[]
    for i in 1:iter
        (length(c42[i]) > 0) && (append!!(c_temp, c42[i]))
    end
    c42 = copy(c_temp)
    
    if length(c1) > 1 && length(c1) < 50_000
        c1 = clique_domination(c1, numthreads)
    end
    
    if length(c2) > 1 && length(c2) < 50_000
        c2 = clique_domination(c2, numthreads)
    end
    
    if length(c3) > 1 && length(c3) < 50_000
        c3 = clique_domination(c3, numthreads)
    end

    if length(c4) > 1 && length(c4) < 50_000
        c4 = clique_domination(c4, numthreads)
    end
    
    return c1, c2, c3, c4, c12, c22, c32, c42
end
