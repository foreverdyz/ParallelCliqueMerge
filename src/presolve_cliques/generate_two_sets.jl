#generate_two_sets.jl

include("build_cg.jl")
include("clique_extend.jl")
include("cliques_dominate.jl")
include("detect_cliques.jl")
include("remap_cons.jl")

function generate_two_set(org_to_bin, bin_to_org, set_pack, set_pack_new, main_cliques, clique_set, var_ub, var_lb, var_type, numthreads)
    cg = build_cg_distributed(set_pack, length(org_to_bin), numthreads)
    cg = update_cg_distributed(set_pack_new, length(org_to_bin), cg, numthreads)
    cg = update_cg_distributed(main_cliques, length(org_to_bin), cg, numthreads)
    cg = update_cg_distributed(clique_set, length(org_to_bin), cg, numthreads)

    c1, c12 = clq_extend_parallel_new(set_pack, cg, n)
    c2, c22 = clq_extend_parallel_new(main_cliques, cg, n)
    c3, c32 = clq_extend_parallel_new(set_pack_new, cg, n)
    c4, c42 = clq_extend_parallel_new(clique_set, cg, n)

    if length(c1) < 100000
        (length(c1) > 1) && (c1 = dominate_cliques_parallel(c1, numthreads))
    end
    if length(c2) < 100000
        (length(c2) > 1) && (c2 = dominate_cliques_parallel(c2, numthreads))
    end
    if length(c3) < 100000
        (length(c3) > 1) && (c3 = dominate_cliques_parallel(c3, numthreads))
    end
    #=
    #since for all cliques from lazy sets and from the no-mian cliques from knapsacks, we will use them as user cut,
    #which will not be added to the original model, but just used to cut off branches
    #So it is unnecessary to chech domination of them
    if length(c4) < 100000
        (length(c4) > 1) && (c4 = dominate_cliques_parallel(c4, numthreads))
    end
    if length(c12) < 100000
        (length(c12) > 1) && (c1 = dominate_cliques_parallel_test(c12, numthreads))
    end
    if length(c22) < 100000
        (length(c22) > 1) && (c2 = dominate_cliques_parallel_test(c22, numthreads))
    end
    if length(c32) < 100000
        (length(c32) > 1) && (c3 = dominate_cliques_parallel_test(c32, numthreads))
    end
    if length(c42) < 100000
        (length(c42) > 1) && (c4 = dominate_cliques_parallel_test(c42, numthreads))
    end
    =#
    
    (length(c1) > 0) && (c1 = remap_cons(c1, bin_to_org, var_lb));
    (length(c2) > 0) && (c2 = remap_cons(c2, bin_to_org, var_lb));
    (length(c3) > 0) && (c3 = remap_cons(c3, bin_to_org, var_lb));
    (length(c4) > 0) && (c4 = remap_cons(c4, bin_to_org, var_lb));
    (length(c12) > 0) && (c12 = remap_cons(c12, bin_to_org, var_lb));
    (length(c22) > 0) && (c22 = remap_cons(c22, bin_to_org, var_lb));
    (length(c32) > 0) && (c32 = remap_cons(c32, bin_to_org, var_lb));
    (length(c42) > 0) && (c42 = remap_cons(c42, bin_to_org, var_lb));
    
    return c1, c2, c3, c4, c12, c22, c32, c42
end