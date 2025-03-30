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
        t = @timed c1, c2, c3, c4, c12, c22, c32, c42 = generate_sets(
            set_pack, set_pack_new, knapsack_set, org_to_bin,
            threadnum, 5_000, 1_000, 25_000_000, 60, 1_250_000, 60
        );

        c1_nnz = 0
        if length(c1) > 0 
            c1_nnz = sum(length(c1[i]) for i in 1:length(c1))
        end
        c2_nnz = 0
        if length(c2) > 0 
            c2_nnz = sum(length(c2[i]) for i in 1:length(c2))
        end
        c3_nnz = 0
        if length(c3) > 0 
            c3_nnz = sum(length(c3[i]) for i in 1:length(c3))
        end
        c4_nnz = 0
        if length(c4) > 0 
            c4_nnz = sum(length(c4[i]) for i in 1:length(c4))
        end
        c12_nnz = 0
        if length(c12) > 0 
            c12_nnz = sum(length(c12[i]) for i in 1:length(c12))
        end
        c22_nnz = 0
        if length(c22) > 0 
            c22_nnz = sum(length(c22[i]) for i in 1:length(c22))
        end
        c32_nnz = 0
        if length(c32) > 0 
            c32_nnz = sum(length(c32[i]) for i in 1:length(c32))
        end
        c42_nnz = 0
        if length(c42) > 0 
            c42_nnz = sum(length(c42[i]) for i in 1:length(c42))
        end
        model  = rebuild_model(
            con_matrix, con_ub, con_lb, var_ub, var_lb, var_type, 
            obj_coef, c1, c2, c3, c4, c12, c22, c32, c42, c1_nnz, c4_nnz, 
            is_min, bin_to_org, filename
        );
        write_to_file(model, "presolved_data_"*string(threadnum)*"/"*filename);
        push!(A, (filename[1:end-4], round(t.time - t.gctime, digits = 4), nnz, setpack_nnz, c1_nnz, c2_nnz, c3_nnz, c4_nnz, c12_nnz, c22_nnz, c32_nnz, c42_nnz))
	println(filename);
    end
    writedlm("presolve_res_"*string(threadnum)*".csv", A, ',');
end
