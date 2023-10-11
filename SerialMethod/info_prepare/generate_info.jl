#generate_info.jl
using BangBang

include("model_info.jl")
include("singleton_and_bounds.jl")
include("binary_list.jl")
include("find_set_pack.jl")
include("constraint_process.jl")

function generate_info(filename::String)
    #get info from model
    con_matrix, con_ub, con_lb, var_ub, var_lb, var_type = model_info(filename)
    #do a simple one round paritial scan: check singleton, and do simple bounds strengthen
    con_matrix, con_ub, con_lb, var_ub, var_lb, var_type = singleton_and_bounds(con_matrix, con_ub, con_lb, var_ub, var_lb, var_type)
    #build two dictionaries
    org_to_bin, bin_to_org = binary_list(var_type)
    #find set packing constraints from the original model
    con_matrix, con_ub, con_lb, set_pack = find_set_pack(con_matrix, con_ub, con_lb, var_type, org_to_bin)
    #process all constraints
    id, set_pack_new, knapsack_set, fix_set = process_all_cons(con_matrix, con_ub, con_lb, var_ub, var_lb, var_type, org_to_bin, bin_to_org)
    if id == false
        println("The problem is infeasible.")
        return false
    end
    #combine two set pack
    (length(set_pack_new) > 0) && (append!!(set_pack, set_pack_new))
    return org_to_bin, bin_to_org, set_pack, knapsack_set, fix_set, con_matrix, con_ub, con_lb, var_ub, var_lb, var_type
end