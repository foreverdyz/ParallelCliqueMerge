#reduced_problem_solve.jl

using Gurobi
using SparseArrays
using Base.Threads
using JuMP

include("info_prepare/generate_info.jl")
include("presolve_cliques/generate_two_sets.jl")
include("rebuild_model.jl")

function reduced_problem_solve(filename::String, solver_threads::Int64)
    #collect necessary info
    org_to_bin, bin_to_org, set_pack, set_pack_new, knapsack_set, fix_set, con_matrix, con_ub, con_lb, var_ub, var_lb, var_type, obj_coef, is_min = generate_info(filename);
    #detect cliques from knapsack constraints
    main_cliques, clique_set = detect_cliques(knapsack_set);
    l1 = length(set_pack)
    l2 = length(set_pack_new)
    l3 = length(knapsack_set)
    @show l1
    @show l2
    @show l3
    #there are too few cliques
    (l1 + l2 + l3 < 3) && (return false)
    
    s1 = 0
    for k in set_pack
        s1 += length(k)
    end
    for k in set_pack_new
        s1 += length(k)
    end
    for k in knapsack_set
        s1 += length(k)
    end

    if s1 > 1_000_000
        println("Give up using cliques methods.")
        return false
    end
    
    println("For the orginal model:")
    model = rebuild_model_org(con_matrix, con_ub, con_lb, var_ub, var_lb, var_type, obj_coef, set_pack, bin_to_org, is_min)
    set_optimizer(model, Gurobi.Optimizer);
    set_silent(model)
    set_attribute(model, "Threads", solver_threads);
    #set_attribute(model, "Presolve", 0);
    set_time_limit_sec(model, 3600.0);
    optimize!(model)
    if has_values(model)
        println("obj is ", objective_value(model))
        println("solve time: ", solve_time(model))
        println("node count: ", node_count(model))
    else
        @show solution_summary(model)
    end
    
    #presolve the problem
    @time c1, c2, c3, c4, c12, c22, c32, c42 = generate_two_set(org_to_bin, bin_to_org, set_pack, set_pack_new, main_cliques, clique_set, var_ub, var_lb, var_type, 24);

    println("For the reduced model:")
    model_rebuild = rebuild_model(con_matrix, con_ub, con_lb, var_ub, var_lb, var_type, obj_coef, c1, c2, c3, c4, c12, c22, c32, c42, is_min, l2, l3);
    set_optimizer(model_rebuild, Gurobi.Optimizer);
    set_silent(model_rebuild)
    set_attribute(model_rebuild, "Threads", solver_threads);
    set_time_limit_sec(model_rebuild, 3600.0);
    optimize!(model_rebuild)
    if has_values(model_rebuild)
        println("obj is ", objective_value(model_rebuild))
        println("solve time: ", solve_time(model_rebuild))
        println("node count: ", node_count(model_rebuild))
    else
        @show solution_summary(model_rebuild)
    end
    
    return true
end