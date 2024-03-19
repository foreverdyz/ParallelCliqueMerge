#presolve_runtime.jl

include("info_prepare/generate_info.jl")
include("clq_merge/generate_sets.jl")
using Base.Threads

"""
we set: 
max_tolerance = 0.00001,
max_length = 3000,
max_nonzero = 1_000_000, 
max_domian = 100_000
"""
function test_runtime(filename::String, numthreads::Int, max_tolerance::Float64, max_length::Int64, max_nonzero::Int64, max_domian::Int64)
    con_matrix, con_set, con_coef, 
    con_ub, con_lb, 
    var_ub, var_lb, var_type, 
    set_pack, set_pack_new, knapsack_set, 
    org_to_bin, bin_to_org, 
    obj_coef, is_min = generate_info("pre_compile_model.mps", max_tolerance);
    c1, c2, c3, c4, c12, c22, c32, c42 = generate_sets(
        set_pack, set_pack_new, knapsack_set, 
        org_to_bin, numthreads
        max_length, max_nonzero, max_domian
    );
    
    con_matrix, con_set, con_coef, 
    con_ub, con_lb, 
    var_ub, var_lb, var_type, 
    set_pack, set_pack_new, knapsack_set, 
    org_to_bin, bin_to_org, 
    obj_coef, is_min = generate_info(filename, max_tolerance);
    t1 = time()
    @time c1, c2, c3, c4, c12, c22, c32, c42 = generate_sets(
        set_pack, set_pack_new, knapsack_set, 
        org_to_bin, numthreads,
        max_length, max_nonzero, max_domian
    );
    t2 = time()
    return filename[9:end-4], numthreads, t2 - t1
end