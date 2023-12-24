#presolve_runtimes.jl

include("info_prepare/generate_info.jl")
include("presolve_cliques/generate_two_sets.jl")
using Random
using Base.Threads
begin
    @info "You are using " * string(nthreads()) * " threads."
    @info "Use \"julia --threads number_threads\" to restart julia and change the number of threads."
end

function clq_merge_runtime(filename::String)
    org_to_bin, bin_to_org, set_pack, set_pack_new, knapsack_set, fix_set, con_matrix, con_ub, con_lb, var_ub, var_lb, var_type, obj_coef, is_min = generate_info(filename)
    println(filename)
    println("Serial Method")
    @time _clq_merge_runtime_serial(org_to_bin, bin_to_org, set_pack, set_pack_new, knapsack_set, var_lb)
    println(string(nthreads())*"Threads")
    @time _clq_merge_runtime_parallel(org_to_bin, bin_to_org, set_pack, set_pack_new, knapsack_set, var_lb, nthreads())
end


function _clq_merge_runtime_serial(org_to_bin, bin_to_org, set_pack, set_pack_new, knapsack_set, var_lb)
    main_cliques, clique_set = detect_cliques(knapsack_set)

    #here we shuffle these set to guarantee both serial and parallel methods are processing same lists (since order may impact the runtime)
    Random.seed!(1234)
    shuffle!(set_pack)
    shuffle!(set_pack_new)
    shuffle!(main_cliques)
    shuffle!(clique_set)

    cg = build_cg_serial(set_pack, length(org_to_bin))
    cg = update_cg_serial(set_pack_new, length(org_to_bin), cg)
    cg = update_cg_serial(main_cliques, length(org_to_bin), cg)
    cg = update_cg_serial(clique_set, length(org_to_bin), cg)

    c1, c12 = cliques_extend_serial(set_pack, cg)
    c2, c22 = cliques_extend_serial(main_cliques, cg)
    c3, c32 = cliques_extend_serial(set_pack_new, cg)
    c4, c42 = cliques_extend_serial(clique_set, cg)

    (length(c1) > 1) && (c1 = dominate_cliques(c1))
    (length(c2) > 1) && (c2 = dominate_cliques(c2))
    (length(c3) > 1) && (c3 = dominate_cliques(c3))
end

function _clq_merge_runtime_parallel(org_to_bin, bin_to_org, set_pack, set_pack_new, knapsack_set, var_lb, n)
    main_cliques, clique_set = detect_cliques_parallel(knapsack_set, n)

    Random.seed!(1234)
    shuffle!(set_pack)
    shuffle!(set_pack_new)
    shuffle!(main_cliques)
    shuffle!(clique_set)

    cg = build_cg_distributed(set_pack, length(org_to_bin), n)
    cg = update_cg_distributed(set_pack_new, length(org_to_bin), cg, n)
    cg = update_cg_distributed(main_cliques, length(org_to_bin), cg, n)
    cg = update_cg_distributed(clique_set, length(org_to_bin), cg, n)

    c1, c12 = clq_extend_parallel_new(set_pack, cg, n)
    c2, c22 = clq_extend_parallel_new(main_cliques, cg, n)
    c3, c32 = clq_extend_parallel_new(set_pack_new, cg, n)
    c4, c42 = clq_extend_parallel_new(clique_set, cg, n)

    (length(c1) > 1) && (c1 = dominate_cliques_parallel(c1, n))
    (length(c2) > 1) && (c2 = dominate_cliques_parallel(c2, n))
    (length(c3) > 1) && (c3 = dominate_cliques_parallel(c3, n))
end
