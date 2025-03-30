#rebuild_model.jl

using BangBang
using JuMP


function rebuild_model(
        con_matrix::AbstractSparseArray, con_ub::AbstractVector, con_lb::AbstractVector,
        var_ub::AbstractVector, var_lb::AbstractVector, var_type::AbstractVector, 
        obj_coef::AbstractVector, 
        c1::Vector{Vector{Int64}}, c2::Vector{Vector{Int64}}, 
        c3::Vector{Vector{Int64}}, c4::Vector{Vector{Int64}},
        c12::Vector{Vector{Int64}}, c22::Vector{Vector{Int64}},
        c32::Vector{Vector{Int64}}, c42::Vector{Vector{Int64}},
        is_min::Bool, bin_to_org::Dict{Int64, Int64}, filename::String,
        nnz::Int64, c2_nnz::Int64, c3_nnz::Int64, c4_nnz::Int64,
    )
    con_num = length(con_ub)
    con_type = zeros(con_num)
    var_num = length(var_ub)
    #denote constraints' types:
    #if lb = -inf, con <= ub
    #if ub = inf, con >= lb
    #otherwise, con == ub or lb
    for index in 1:con_num
        if con_lb[index] == -Inf
            con_type[index] = -1
        elseif con_ub[index] == Inf
            con_type[index] = 1
        end
    end
    for index in 1:con_num
        if con_lb[index] == -Inf
            if con_ub[index] == Inf
                con_type[index] = 2
            end
        end
    end
    
    model_rebuild = Model();
    #build variables
    @variable(model_rebuild, var_ub[i] >= x[i in 1:var_num] >= var_lb[i])
    for i in 1:var_num
        (var_type[i] == 2) && (set_binary(x[i]))
        (var_type[i] == 1) && (set_integer(x[i]))
    end

    #build constraints based on constraints' types
    for j in 1:con_num
        if con_type[j] == 1
            @constraint(model_rebuild, sum(con_matrix[j, i] * x[i] for i in findnz(con_matrix[j, :])[1]) >= con_lb[j])
        elseif con_type[j] == -1
            @constraint(model_rebuild, sum(con_matrix[j, i] * x[i] for i in findnz(con_matrix[j, :])[1]) <= con_ub[j])
        elseif  con_type[j] == 0
            @constraint(model_rebuild, sum(con_matrix[j, i] * x[i] for i in findnz(con_matrix[j, :])[1]) == con_lb[j])
        end
    end
    
    if length(c1) > 0
        c1 = _remap_cons(c1, bin_to_org, var_lb)
        @constraint(
            model_rebuild,
            con1[j in 1:length(c1)],
            sum(c1[j][i] * x[i] for i in findnz(c1[j])[1]) <= 1 + sum(min.(findnz(c1[j])[2], 0))
        )
    end
    
    if length(c2) > 0
        c2 = _remap_cons(c2, bin_to_org, var_lb)
        if c2_nnz/nnz < 1
            @constraint(
                model_rebuild,
                con2[j in 1:length(c2)],
                sum(c2[j][i] * x[i] for i in findnz(c2[j])[1]) <= 1 + sum(min.(findnz(c2[j])[2], 0))
            )
        else
            @constraint(
                model_rebuild,
                con2[j in 1:length(c2)],
                sum(c2[j][i] * x[i] for i in findnz(c2[j])[1]) <= 1 + sum(min.(findnz(c2[j])[2], 0)),
                base_name = "lazy_c2"
            )
        end
    end
    
    if length(c3) > 0
        c3 = _remap_cons(c3, bin_to_org, var_lb)
        if c3_nnz/nnz < 1
            @constraint(
                model_rebuild,
                con3[j in 1:length(c3)],
                sum(c3[j][i] * x[i] for i in findnz(c3[j])[1]) <= 1 + sum(min.(findnz(c3[j])[2], 0))
            )
        else
            @constraint(
                model_rebuild,
                con3[j in 1:length(c3)],
                sum(c3[j][i] * x[i] for i in findnz(c3[j])[1]) <= 1 + sum(min.(findnz(c3[j])[2], 0)),
                base_name = "lazy_c3"
            )
        end
    end

    if length(c4) > 0
        c4 = _remap_cons(c4, bin_to_org, var_lb)
        if (c2_nnz + c4_nnz)/nnz < 1
            @constraint(
                model_rebuild,
                con4[j in 1:length(c4)],
                sum(c4[j][i] * x[i] for i in findnz(c4[j])[1]) <= 1 + sum(min.(findnz(c4[j])[2], 0))
            )
        else
            @constraint(
                model_rebuild,
                con4[j in 1:length(c4)],
                sum(c4[j][i] * x[i] for i in findnz(c4[j])[1]) <= 1 + sum(min.(findnz(c4[j])[2], 0)),
                base_name = "lazy_c4"
            )
        end
    end
    
    if length(c12) > 0
        c12 = _remap_cons(c12, bin_to_org, var_lb)
        @constraint(
            model_rebuild,
            con12[j in 1:length(c12)],
            sum(c12[j][i] * x[i] for i in findnz(c12[j])[1]) <= 1 + sum(min.(findnz(c12[j])[2], 0)),
            base_name = "lazy_c12"
        )
    end
    
    if length(c22) > 0
        c22 = _remap_cons(c22, bin_to_org, var_lb)
        @constraint(
            model_rebuild,
            con22[j in 1:length(c22)],
            sum(c22[j][i] * x[i] for i in findnz(c22[j])[1]) <= 1 + sum(min.(findnz(c22[j])[2], 0)),
            base_name = "lazy_c22"
        )
    end
    
    if length(c32) > 0
        c32 = _remap_cons(c32, bin_to_org, var_lb)
        @constraint(
            model_rebuild,
            con32[j in 1:length(c32)],
            sum(c32[j][i] * x[i] for i in findnz(c32[j])[1]) <= 1 + sum(min.(findnz(c32[j])[2], 0)),
            base_name = "lazy_c32"
        )
    end
    
    if length(c42) > 0
        c42 = _remap_cons(c42, bin_to_org, var_lb)
        @constraint(
            model_rebuild,
            con42[j in 1:length(c42)],
            sum(c42[j][i] * x[i] for i in findnz(c42[j])[1]) <= 1 + sum(min.(findnz(c42[j])[2], 0)),
            base_name = "lazy_c42"
        )
    end
    
    #set objective function
    (is_min) ? (@objective(model_rebuild, Min, sum(obj_coef[i] * x[i] for i in 1:var_num))) : (@objective(model_rebuild, Max, sum(obj_coef[i] * x[i] for i in 1:var_num)))
    
    
    return model_rebuild
end


function _remap_cons(set::Vector{Vector{Int64}}, bin_to_org::Dict{Int64, Int64}, var_lb::AbstractVector)
    m = length(set)
    gen_set = []
    binary_number = length(bin_to_org)
    for j in 1:m
        c = spzeros(length(var_lb))
        for i in set[j] 
            (i <= binary_number) ? (c[bin_to_org[i]] = 1) : (c[bin_to_org[i - binary_number]] = -1)
        end
        push!(gen_set, c)
    end
    return gen_set
end
#=
#rebuild_model.jl

#ENV["GRB_LICENSE_FILE"] = "/usr/local/gurobi/10.0.1/gurobi.lic"
#using Gurobi
using BangBang
using JuMP


function rebuild_model(
        con_matrix::AbstractSparseArray, con_ub::AbstractVector, con_lb::AbstractVector,
        var_ub::AbstractVector, var_lb::AbstractVector, var_type::AbstractVector,
        obj_coef::AbstractVector,
        c1::Vector{Vector{Int64}}, c2::Vector{Vector{Int64}},
        c3::Vector{Vector{Int64}}, c4::Vector{Vector{Int64}},
        c12::Vector{Vector{Int64}}, c22::Vector{Vector{Int64}},
        c32::Vector{Vector{Int64}}, c42::Vector{Vector{Int64}},
        nnz::Int64,
        is_min::Bool, bin_to_org::Dict{Int64, Int64}, filename::String
    )
    alpha = 0.3;
    con_num = length(con_ub)
    con_type = zeros(con_num)
    var_num = length(var_ub)
    #denote constraints' types:
    #if lb = -inf, con <= ub
    #if ub = inf, con >= lb
    #otherwise, con == ub or lb
    for index in 1:con_num
        if con_lb[index] == -Inf
            con_type[index] = -1
        elseif con_ub[index] == Inf
            con_type[index] = 1
        end
    end
    for index in 1:con_num
        if con_lb[index] == -Inf
            if con_ub[index] == Inf
                con_type[index] = 2
            end
        end
    end

    model_rebuild = Model();
    #build variables
    @variable(model_rebuild, var_ub[i] >= x[i in 1:var_num] >= var_lb[i])
    for i in 1:var_num
        (var_type[i] == 2) && (set_binary(x[i]))
        (var_type[i] == 1) && (set_integer(x[i]))
    end

    #build constraints based on constraints' types
    for j in 1:con_num
        if con_type[j] == 1
            @constraint(model_rebuild, sum(con_matrix[j, i] * x[i] for i in findnz(con_matrix[j, :])[1]) >= con_lb[j])
        elseif con_type[j] == -1
            @constraint(model_rebuild, sum(con_matrix[j, i] * x[i] for i in findnz(con_matrix[j, :])[1]) <= con_ub[j])
        elseif  con_type[j] == 0
            @constraint(model_rebuild, sum(con_matrix[j, i] * x[i] for i in findnz(con_matrix[j, :])[1]) == con_lb[j])
        end
    end

    if length(c1) > 0
        c1 = _remap_cons(c1, bin_to_org, var_lb)
        @constraint(
            model_rebuild,
            con1[j in 1:length(c1)],
            sum(c1[j][i] * x[i] for i in findnz(c1[j])[1]) <= 1 + sum(min.(findnz(c1[j])[2], 0)),
            base_name = "lazy_c1"
        )
    end

    c2_nnz = 0
    if length(c2) > 0
        c2 = _remap_cons(c2, bin_to_org, var_lb)
        c2_nnz = sum(length(c2[i]) for i in 1:length(c2))
        @constraint(
            model_rebuild,
            con2[j in 1:length(c2)],
            sum(c2[j][i] * x[i] for i in findnz(c2[j])[1]) <= 1 + sum(min.(findnz(c2[j])[2], 0)),
            base_name = "lazy_c2"
        )
        if c2_nnz >= alpha*nnz
            for i in 1:length(c2)
                MOI.set(model_rebuild, Gurobi.ConstraintAttribute("Lazy"), con2[i], -1)
            end
        else
            c2_nnz = 0
        end
    end

    c3_nnz = 0
    if length(c3) > 0
        c3 = _remap_cons(c3, bin_to_org, var_lb)
        c3_nnz = sum(length(c3[i]) for i in 1:length(c3))
        @constraint(
            model_rebuild,
            con3[j in 1:length(c3)],
            sum(c3[j][i] * x[i] for i in findnz(c3[j])[1]) <= 1 + sum(min.(findnz(c3[j])[2], 0)),
            base_name = "lazy_c3"
        )
        if c2_nnz + c3_nnz >= alpha*nnz
            for i in 1:length(c3)
                MOI.set(model_rebuild, Gurobi.ConstraintAttribute("Lazy"), con3[i], -1)
            end
        else
            c3_nnz = 0
        end
    end

    c4_nnz = 0
    if length(c4) > 0
        c4 = _remap_cons(c4, bin_to_org, var_lb)
        c4_nnz = sum(length(c4[i]) for i in 1:length(c4))
        @constraint(
            model_rebuild,
            con4[j in 1:length(c4)],
            sum(c4[j][i] * x[i] for i in findnz(c4[j])[1]) <= 1 + sum(min.(findnz(c4[j])[2], 0)),
            base_name = "lazy_c4"
        )
        if c2_nnz + c3_nnz + c4_nnz >= alpha*nnz
            for i in 1:length(c4)
                MOI.set(model_rebuild, Gurobi.ConstraintAttribute("Lazy"), con4[i], -1)
            end
        end
    end

    if length(c12) > 0
        c12 = _remap_cons(c12, bin_to_org, var_lb)
        @constraint(
            model_rebuild,
            con12[j in 1:length(c12)],
            sum(c12[j][i] * x[i] for i in findnz(c12[j])[1]) <= 1 + sum(min.(findnz(c12[j])[2], 0)),
            base_name = "lazy_c12"
        )
        for i in 1:length(c12)
            MOI.set(model_rebuild, Gurobi.ConstraintAttribute("Lazy"), con12[i], -1)
        end
    end

    if length(c22) > 0
        c22 = _remap_cons(c22, bin_to_org, var_lb)
        @constraint(
            model_rebuild,
            con22[j in 1:length(c22)],
            sum(c22[j][i] * x[i] for i in findnz(c22[j])[1]) <= 1 + sum(min.(findnz(c22[j])[2], 0)),
            base_name = "lazy_c22"
        )
        for i in 1:length(c22)
            MOI.set(model_rebuild, Gurobi.ConstraintAttribute("Lazy"), con22[i], -1)
        end
    end

    if length(c32) > 0
        c32 = _remap_cons(c32, bin_to_org, var_lb)
        @constraint(
            model_rebuild,
            con32[j in 1:length(c32)],
            sum(c32[j][i] * x[i] for i in findnz(c32[j])[1]) <= 1 + sum(min.(findnz(c32[j])[2], 0)),
            base_name = "lazy_c32"
        )
        for i in 1:length(c32)
            MOI.set(model_rebuild, Gurobi.ConstraintAttribute("Lazy"), con32[i], -1)
        end
    end

    if length(c42) > 0
        c42 = _remap_cons(c42, bin_to_org, var_lb)
        @constraint(
            model_rebuild,
            con42[j in 1:length(c42)],
            sum(c42[j][i] * x[i] for i in findnz(c42[j])[1]) <= 1 + sum(min.(findnz(c42[j])[2], 0)),
            base_name = "lazy_c42"
        )
        for i in 1:length(c42)
            MOI.set(model_rebuild, Gurobi.ConstraintAttribute("Lazy"), con42[i], -1)
        end
    end

    #set objective function
    (is_min) ? (@objective(model_rebuild, Min, sum(obj_coef[i] * x[i] for i in 1:var_num))) : (@objective(model_rebuild, Max, sum(obj_coef[i] * x[i] for i in 1:var_num)))

    return model_rebuild
end


function _remap_cons(set::Vector{Vector{Int64}}, bin_to_org::Dict{Int64, Int64}, var_lb::AbstractVector)
    m = length(set)
    gen_set = []
    binary_number = length(bin_to_org)
    for j in 1:m
        c = spzeros(length(var_lb))
        for i in set[j]
            (i <= binary_number) ? (c[bin_to_org[i]] = 1) : (c[bin_to_org[i - binary_number]] = -1)
        end
        push!(gen_set, c)
    end
    return gen_set
end
=#
