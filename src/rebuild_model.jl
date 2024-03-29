#rebuild_model.jl

using Gurobi
using BangBang
using JuMP

"""
    rebuild_model_0(): build model only with org set pack
    rebuild_model_1(): build model with all sets we generated
    rebuild_model_2(): build direct_model only with org set pack
    rebuild_model_3(): build direct_model with all sets we generated
"""
function rebuild_model_0(
        con_matrix::AbstractSparseArray, con_ub::AbstractVector, con_lb::AbstractVector,
        var_ub::AbstractVector, var_lb::AbstractVector, var_type::AbstractVector, 
        obj_coef::AbstractVector, 
        c1::Vector{Vector{Int64}}, 
        is_min::Bool, bin_to_org::Dict{Int64, Int64}
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
    
    #check we get correct info of model
    #model_rebuild = direct_model(Gurobi.Optimizer())#Gurobi.Optimizer)
    model_rebuild = Model();
    #build variables
    @variable(model_rebuild, var_ub[i] >= x[i in 1:var_num] >= var_lb[i])

    #set binary and integer variables
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
    
    #set objective function
    (is_min) ? (@objective(model_rebuild, Min, sum(obj_coef[i] * x[i] for i in 1:var_num))) : (@objective(model_rebuild, Max, sum(obj_coef[i] * x[i] for i in 1:var_num)))
    
    
    return model_rebuild
end

function rebuild_model_1(
        con_matrix::AbstractSparseArray, con_ub::AbstractVector, con_lb::AbstractVector,
        var_ub::AbstractVector, var_lb::AbstractVector, var_type::AbstractVector, 
        obj_coef::AbstractVector, 
        c1::Vector{Vector{Int64}}, c2::Vector{Vector{Int64}}, 
        c3::Vector{Vector{Int64}}, c4::AbstractVector,
        c12::AbstractVector, c22::AbstractVector,
        c32::AbstractVector, c42::AbstractVector, 
        is_min::Bool, bin_to_org::Dict{Int64, Int64},
        c2_final::Bool, c3_final::Bool
    )
    #combine cliques from c12 to one set
    iter = length(c12)
    c_temp = Vector{Int64}[]
    for i in 1:iter
        (length(c12[i]) > 0) && (append!!(c_temp, c12[i]))
    end
    c12 = copy(c_temp)
    iter = length(c22)
    c_temp = Vector{Int64}[]
    for i in 1:iter
        (length(c22[i]) > 0) && (append!!(c_temp, c22[i]))
    end
    c22 = copy(c_temp)
    iter = length(c32)
    c_temp = Vector{Int64}[]
    for i in 1:iter
        (length(c32[i]) > 0) && (append!!(c_temp, c32[i]))
    end
    c32 = copy(c_temp)

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
    
    #check we get correct info of model
    #model_rebuild = direct_model(Gurobi.Optimizer())#Gurobi.Optimizer)
    model_rebuild = Model();
    #build variables
    @variable(model_rebuild, var_ub[i] >= x[i in 1:var_num] >= var_lb[i])

    #set binary and integer variables
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
        #for i in 1:length(c1)
        #    MOI.set(model_rebuild, Gurobi.ConstraintAttribute("Lazy"), con1[i], 3)
        #end
    end
    
    if length(c2) > 0
        c2 = _remap_cons(c2, bin_to_org, var_lb)
        @constraint(
            model_rebuild,
            con2[j in 1:length(c2)],
            sum(c2[j][i] * x[i] for i in findnz(c2[j])[1]) <= 1 + sum(min.(findnz(c2[j])[2], 0))
        )
        if c2_final
            for i in 1:length(c2)
                MOI.set(model_rebuild, Gurobi.ConstraintAttribute("Lazy"), con2[i], -1)
            end
        end
    end
    
    if length(c3) > 0
        c3 = _remap_cons(c3, bin_to_org, var_lb)
        @constraint(
            model_rebuild,
            con3[j in 1:length(c3)],
            sum(c3[j][i] * x[i] for i in findnz(c3[j])[1]) <= 1 + sum(min.(findnz(c3[j])[2], 0))
        )
        if c3_final
            for i in 1:length(c3)
                MOI.set(model_rebuild, Gurobi.ConstraintAttribute("Lazy"), con3[i], -1)
            end
        end
    end
    
    if length(c4) > 0
        c4 = _remap_cons(c4, bin_to_org, var_lb)
        @constraint(
            model_rebuild,
            con4[j in 1:length(c4)],
            sum(c4[j][i] * x[i] for i in findnz(c4[j])[1]) <= 1 + sum(min.(findnz(c4[j])[2], 0))
        )
        for i in 1:length(c4)
            MOI.set(model_rebuild, Gurobi.ConstraintAttribute("Lazy"), con4[i], -1)
        end
       
    end
    
    if length(c12) > 0
        c12 = _remap_cons(c12, bin_to_org, var_lb)
        @constraint(
            model_rebuild,
            con12[j in 1:length(c12)],
            sum(c12[j][i] * x[i] for i in findnz(c12[j])[1]) <= 1 + sum(min.(findnz(c12[j])[2], 0))
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
            sum(c22[j][i] * x[i] for i in findnz(c22[j])[1]) <= 1 + sum(min.(findnz(c22[j])[2], 0))
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
            sum(c32[j][i] * x[i] for i in findnz(c32[j])[1]) <= 1 + sum(min.(findnz(c32[j])[2], 0))
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
            sum(c42[j][i] * x[i] for i in findnz(c42[j])[1]) <= 1 + sum(min.(findnz(c42[j])[2], 0))
        )
       
        for i in 1:length(c42)
            MOI.set(model_rebuild, Gurobi.ConstraintAttribute("Lazy"), con42[i], -1)
        end 
    end
    
    
    #set objective function
    (is_min) ? (@objective(model_rebuild, Min, sum(obj_coef[i] * x[i] for i in 1:var_num))) : (@objective(model_rebuild, Max, sum(obj_coef[i] * x[i] for i in 1:var_num)))
    
    
    return model_rebuild
end

function rebuild_model_2(
        con_matrix::AbstractSparseArray, con_ub::AbstractVector, con_lb::AbstractVector,
        var_ub::AbstractVector, var_lb::AbstractVector, var_type::AbstractVector, 
        obj_coef::AbstractVector, 
        c1::Vector{Vector{Int64}}, 
        is_min::Bool, bin_to_org::Dict{Int64, Int64}
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
    
    #check we get correct info of model
    model_rebuild = direct_model(Gurobi.Optimizer())#Gurobi.Optimizer)
    #model_rebuild = Model();
    #build variables
    @variable(model_rebuild, var_ub[i] >= x[i in 1:var_num] >= var_lb[i])

    #set binary and integer variables
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
    
    #set objective function
    (is_min) ? (@objective(model_rebuild, Min, sum(obj_coef[i] * x[i] for i in 1:var_num))) : (@objective(model_rebuild, Max, sum(obj_coef[i] * x[i] for i in 1:var_num)))
    
    
    return model_rebuild
end

function rebuild_model_3(
        con_matrix::AbstractSparseArray, con_ub::AbstractVector, con_lb::AbstractVector,
        var_ub::AbstractVector, var_lb::AbstractVector, var_type::AbstractVector, 
        obj_coef::AbstractVector, 
        c1::Vector{Vector{Int64}}, c2::Vector{Vector{Int64}}, 
        c3::Vector{Vector{Int64}}, c4::AbstractVector,
        c12::AbstractVector, c22::AbstractVector,
        c32::AbstractVector, c42::AbstractVector, 
        is_min::Bool, bin_to_org::Dict{Int64, Int64},
        c2_final::Bool, c3_final::Bool
    )
    #combine cliques from c12 to one set
    iter = length(c12)
    c_temp = Vector{Int64}[]
    for i in 1:iter
        (length(c12[i]) > 0) && (append!!(c_temp, c12[i]))
    end
    c12 = copy(c_temp)
    iter = length(c22)
    c_temp = Vector{Int64}[]
    for i in 1:iter
        (length(c22[i]) > 0) && (append!!(c_temp, c22[i]))
    end
    c22 = copy(c_temp)
    iter = length(c32)
    c_temp = Vector{Int64}[]
    for i in 1:iter
        (length(c32[i]) > 0) && (append!!(c_temp, c32[i]))
    end
    c32 = copy(c_temp)
    
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
    
    #check we get correct info of model
    model_rebuild = direct_model(Gurobi.Optimizer())#Gurobi.Optimizer)
    #model_rebuild = Model();
    #build variables
    @variable(model_rebuild, var_ub[i] >= x[i in 1:var_num] >= var_lb[i])

    #set binary and integer variables
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
        #for i in 1:length(c1)
        #    MOI.set(model_rebuild, Gurobi.ConstraintAttribute("Lazy"), con1[i], 3)
        #end
    end
    
    if length(c2) > 0
        c2 = _remap_cons(c2, bin_to_org, var_lb)
        @constraint(
            model_rebuild,
            con2[j in 1:length(c2)],
            sum(c2[j][i] * x[i] for i in findnz(c2[j])[1]) <= 1 + sum(min.(findnz(c2[j])[2], 0))
        )
        if c2_final
            for i in 1:length(c2)
                MOI.set(model_rebuild, Gurobi.ConstraintAttribute("Lazy"), con2[i], -1)
            end
        end
    end
    
    if length(c3) > 0
        c3 = _remap_cons(c3, bin_to_org, var_lb)
        @constraint(
            model_rebuild,
            con3[j in 1:length(c3)],
            sum(c3[j][i] * x[i] for i in findnz(c3[j])[1]) <= 1 + sum(min.(findnz(c3[j])[2], 0))
        )
        if c3_final
            for i in 1:length(c3)
                MOI.set(model_rebuild, Gurobi.ConstraintAttribute("Lazy"), con3[i], -1)
            end
        end
    end
    
    if length(c4) > 0
        c4 = _remap_cons(c4, bin_to_org, var_lb)
        @constraint(
            model_rebuild,
            con4[j in 1:length(c4)],
            sum(c4[j][i] * x[i] for i in findnz(c4[j])[1]) <= 1 + sum(min.(findnz(c4[j])[2], 0))
        )
        for i in 1:length(c4)
            MOI.set(model_rebuild, Gurobi.ConstraintAttribute("Lazy"), con4[i], -1)
        end
       
    end
    
    if length(c12) > 0
        c12 = _remap_cons(c12, bin_to_org, var_lb)
        @constraint(
            model_rebuild,
            con12[j in 1:length(c12)],
            sum(c12[j][i] * x[i] for i in findnz(c12[j])[1]) <= 1 + sum(min.(findnz(c12[j])[2], 0))
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
            sum(c22[j][i] * x[i] for i in findnz(c22[j])[1]) <= 1 + sum(min.(findnz(c22[j])[2], 0))
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
            sum(c32[j][i] * x[i] for i in findnz(c32[j])[1]) <= 1 + sum(min.(findnz(c32[j])[2], 0))
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
            sum(c42[j][i] * x[i] for i in findnz(c42[j])[1]) <= 1 + sum(min.(findnz(c42[j])[2], 0))
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
