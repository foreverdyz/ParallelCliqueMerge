#generate_info.jl

include("model_info.jl");
include("simple_presolve.jl");
include("find_setpack.jl");
include("process_all_cons.jl");

using BangBang

function generate_info(filename::String)
    con_matrix, con_set, con_coef, con_ub, con_lb, var_ub, var_lb, 
    var_name, var_type, obj_coef, obj_constant, is_min, nnz = model_info(filename);
    con_matrix, con_set, con_coef, con_ub, con_lb, 
    var_ub, var_lb, var_type, org_to_bin, bin_to_org = simple_presolve(
        con_matrix, con_set, con_coef, 
        con_ub, con_lb, var_ub, var_lb, var_type
    );
    
    I, J, set_pack, con_ub, con_lb = find_setpack(con_set, con_coef, con_ub, con_lb, var_type, org_to_bin);
    setpack_nnz = 0
    for s in set_pack
        setpack_nnz += length(s)
    end

    set_pack_new, knapsack_set = process_all_cons(con_set[I], con_coef[I], con_ub[I], con_lb[I], var_ub, var_lb, var_type, org_to_bin);
    
    #J includes constraints with equality and setpacking (<= or >=)
    if !isempty(J)
        append!!(I, J) #add J back
        sort!(I) #keep the original constraints order
    end
    con_matrix, con_ub, con_lb, con_set, con_coef = con_matrix[I, :], con_ub[I], con_lb[I], con_set[I], con_coef[I];
    
    return con_matrix, con_set, con_coef, con_ub, con_lb, var_ub, var_lb, var_type, set_pack, set_pack_new, knapsack_set, org_to_bin, bin_to_org, obj_coef, obj_constant, is_min, nnz, setpack_nnz
end