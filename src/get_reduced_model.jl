#model_runtimes.jl

include("info_prepare/generate_info.jl")
include("clq_merge/generate_sets.jl")
include("model_build.jl")
using BangBang
using Random

"""
    clq_detect_time_limit: time limit for clq extension
    cg_length_limit: size limit of clqs in cg construction
    cg_term_limit: number of nonzero terms limit in one thread!!! in cg construction
    clq_ext_length_limit: size limit of clqs in clq extension
    clq_ext_term_limit: number of nonzero terms limit in one thread!!! in clq extension
    domain_limit: number of constraints limit in clq domain
"""
function get_reduced_model(
    filename::String, numthreads::Int64,
    clq_detect_time_limit::Real, cg_length_limit::Int64, cg_term_limit::Int64,
    clq_ext_length_limit::Int64, clq_ext_term_limit::Int64, domain_limit::Int64
    )
    #since we need to shuffle clqs in parallel computing, we set the random seed for reproducing
    Random.seed!(12)
    #read model info from the .mps file
    con_matrix, con_set, con_coef, 
    con_ub, con_lb, 
    var_ub, var_lb, var_type, 
    set_pack, set_pack_new, knapsack_set, 
    org_to_bin, bin_to_org, 
    obj_coef, is_min = generate_info(filename);
    #here Gurobi will do the same thing (clq mering in Gurobi)
    if length(set_pack_new) + length(knapsack_set) < 1
        model = read_from_file(filename)
        return model
    end
    
    #report the runtime of our cg procedure
    @time c1, c2, c3, c4, c12, c22, c32, c42 = generate_sets(
        set_pack, set_pack_new, knapsack_set, org_to_bin, numthreads,
        clq_detect_time_limit, cg_length_limit, cg_term_limit,
        clq_ext_length_limit, clq_ext_term_limit, domain_limit
        );
    
    #judge whether clqs from c2 and c3 should be placed as user cuts
    c2_user_cut = false
    if length(c1) < 20_000 && length(c1) >= 120
        if length(c1) > 2000
            (length(c2)/length(c1) >= 0.02) && (c2_final = true)
        else
            (length(c2)/length(c1) >= 1) && (c2_final = true)
        end
    end

    c3_user_cut = false
    if length(c1) < 12_000 && length(c1) >= 1500
        c3_user_cut = true
    end
    
    model = model_build(
        con_matrix, con_ub, con_lb, 
        var_ub, var_lb, var_type, obj_coef, 
        c1, c2, c3, c4, c12, c22, c32, c42, 
        is_min, bin_to_org,
        c2_user_cut, c3_user_cut
    );
    
    return model 
    #!!!!here we cannot use write_to_file() to save the problem as mps file
    #since julia cannot write file with gurobi callback function (user cut/lazy cut)
end