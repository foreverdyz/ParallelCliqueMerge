#build_model_other_solvers.jl

include("rebuild_model_general.jl")
include("info_prepare/generate_info.jl")
include("clq_merge/generate_sets.jl")
using BangBang
using Random

"""
we set: 
max_tolerance = 0.00001,
max_length = 3000,
max_nonzero = 1_000_000, 
max_domian = 100_000

All constraints from generated .mps file starting with "user_lazy_con" should be placed as user cuts.
"""
function model_mps(filename::String, numthreads::Int64, max_tolerance::Float64, max_length::Int64, max_nonzero::Int64, max_domian::Int64)
    Random.seed!(12)
    con_matrix, con_set, con_coef, 
    con_ub, con_lb, 
    var_ub, var_lb, var_type, 
    set_pack, set_pack_new, knapsack_set, 
    org_to_bin, bin_to_org, 
    obj_coef, is_min = generate_info(filename, max_tolerance);
    if length(set_pack_new) + length(knapsack_set) < 1
        model_rebuild = rebuild_model_0(
            con_matrix, con_ub, con_lb, 
            var_ub, var_lb, var_type, obj_coef, 
            set_pack,
            is_min, bin_to_org
         );
        write_to_file(model_rebuild, filename[9:end-4]*".mps")
        return "Completed"
    end
    
    c1, c2, c3, c4, c12, c22, c32, c42 = generate_sets(
        set_pack, set_pack_new, knapsack_set, 
        org_to_bin, numthreads, 
        max_length, max_nonzero, max_domian
    );
    
    #c2_final or c3_final = true mean this set will be placed as user cuts
    c2_final = true
    c3_final = true
    
    if length(c1) > 20_000
        c2_final = false
        c3_final = false
    end
    
    if length(c1) < 50
        c2_final = false
        c3_final = false
    end
    
    (length(c3)/ length(c1) < 0.005) && (c3_final = false)
    
    if length(c1) > 10_000 && length(c1) <= 20_000
        (length(c2)/ length(c1) < 0.02) && (c2_final = false)
    end
    
    if length(c1) > 100 && length(c1) <= 10_000
        (length(c2)/length(c1) < 1) && (c2_final = false)
    end
    
    if length(c2) > 8000
        c2_final = true
    end
    
    if length(c3) > 8000
        c3_final = true
    end
    
    model_rebuild = rebuild_model_1(
        con_matrix, con_ub, con_lb, 
        var_ub, var_lb, var_type, obj_coef, 
        c1, c2, c3, c4, c12, c22, c32, c42, 
        is_min, bin_to_org,
        c2_final, c3_final
    );
    write_to_file(model_rebuild, filename[9:end-4]*".mps")
    return "Completed"
end
