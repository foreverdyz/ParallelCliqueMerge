#generate_info_new.jl

include("model_info.jl");
include("simple_presolve.jl");
include("binary_list.jl");
include("find_setpack.jl");
include("process_all_cons.jl")

using BangBang
"""
    generate_info()
    includes reading info from file, a one-round simple presolve, and detect original set packs, implied set packs, and conflicting knapsacks 
"""
function generate_info(filename::String, max_tolerance::Float64)
    #read from file
    con_matrix, con_set, con_coef, con_ub, con_lb, var_ub, var_lb, var_type, obj_coef, is_min = model_info(filename);
    #one-round simple presolve
    con_matrix, con_set, con_coef, 
    con_ub, con_lb, var_ub, var_lb, var_type = simple_presolve(
        con_matrix, con_set, con_coef, 
        con_ub, con_lb, var_ub, var_lb, 
        var_type, max_tolerance
    );
    #build dictionary for nonbinary to binary and binary to nonbinary
    org_to_bin, bin_to_org = binary_list(var_type);
    #find original set packs
    I, J, set_pack, con_ub, con_lb = find_set_pack(
        con_set, con_coef, 
        con_ub, con_lb,
        var_type, org_to_bin
    );
    #find implied set packs and conflicting knapsacks
    set_pack_new, knapsack_set = process_all_cons_new(
        con_set[I], con_coef[I],
        con_ub[I], con_lb[I],
        var_ub, var_lb,
        var_type, org_to_bin
    );
    #I and J includes all indeces of constraints have not been removed from MIP
    (length(J) > 0) && (append!!(I, J))
    con_matrix, con_ub, con_lb, con_set, con_coef = con_matrix[I, :], con_ub[I], con_lb[I], con_set[I], con_coef[I];
    return con_matrix, con_set, con_coef, con_ub, con_lb, var_ub, var_lb, var_type, set_pack, set_pack_new, knapsack_set, org_to_bin, bin_to_org, obj_coef, is_min
end