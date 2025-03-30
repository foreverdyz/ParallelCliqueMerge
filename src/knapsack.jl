#presolve.jl

include("prepresolves/generate_info.jl");
include("generate_cliques/generate_sets.jl");
include("build_model/rebuild_model.jl");

using DelimitedFiles

list_1 = open("first_round.txt") do f
    readlines(f)
end;

threadnum = parse(Int64, ARGS[1]);

if threadnum > 0
    #collect all res with A
    A = [];
    for i in 1:length(list_1)
        filename = list_1[i];
        con_matrix, con_set, con_coef, con_ub, con_lb, 
        var_ub, var_lb, var_type, set_pack, set_pack_new, knapsack_set, 
        org_to_bin, bin_to_org, obj_coef, obj_constant, is_min, 
        nnz, setpack_nnz = generate_info("mipdata/"*filename);
        #release memory
        con_set, con_coef = [], []
        t = generate_sets(
            set_pack, set_pack_new, knapsack_set, org_to_bin,
            threadnum, 1_000, 1_000, 25_000_000, 60, 1_250_000, 60
        ); #in our own pc (64GB), we set 5_000, 1_000, 100_000_000, 60, 5_000_000, 60
	println(filename)
	println(t)
    end
end
