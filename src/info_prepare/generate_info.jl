#generate_info_new.jl

include("model_info.jl");
include("simple_presolve.jl");
include("binary_list.jl");
include("find_setpack.jl");
include("process_all_cons.jl")

using BangBang

"""
    generate_info() reads model from the .mps file, conducts one round simple presolve, 
    and detects set packing constraints and conflicted knapsacks.
"""
function generate_info(filename::String)
    #read model from .mps file
    con_matrix, con_set, con_coef, con_ub, con_lb, var_ub, var_lb, var_type, obj_coef, is_min = model_info(filename);
    
    #empty and singletons remvoing, and one round varaibles bounds strengthening (domain propagation)
    con_matrix, con_set, con_coef, 
    con_ub, con_lb, var_ub, var_lb, var_type = simple_presolve(
        con_matrix, con_set, con_coef, 
        con_ub, con_lb, var_ub, var_lb, var_type
    );
    
    #build dictionaries for org variables to bin variables and conversely 
    org_to_bin, bin_to_org = binary_list(var_type);

    #detect set packing (cliques) constraints
    I, J, set_pack, con_ub, con_lb = find_set_pack(
        con_set, con_coef, 
        con_ub, con_lb,
        var_type, org_to_bin
    );

    #detect new set packing constraints and conflicted knapsacks (in which we can detect new clqs later)
    set_pack_new, knapsack_set = process_all_cons_new(
        con_set[I], con_coef[I],
        con_ub[I], con_lb[I],
        var_ub, var_lb,
        var_type, org_to_bin
    );

    #I: non set packing constraints, J: set partition constraints (treating them as two constraints, and one of thme is a clq)
    (length(J) > 0) && (append!!(I, J))
    con_matrix, con_ub, con_lb, con_set, con_coef = con_matrix[I, :], con_ub[I], con_lb[I], con_set[I], con_coef[I];

    return con_matrix, con_set, con_coef, con_ub, con_lb, var_ub, var_lb, var_type, set_pack, set_pack_new, knapsack_set, org_to_bin, bin_to_org, obj_coef, is_min
end