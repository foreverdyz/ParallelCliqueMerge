#rebuild_model.jl

function rebuild_model(con_matrix, con_ub, con_lb, var_ub, var_lb, var_type, obj_coef, c1, c2, c3, c4, c12, c22, c32, c42, is_min, l2, l3)
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
    model_rebuild = Model()#Gurobi.Optimizer)

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
        @constraint(
            model_rebuild,
            con1[j in 1:length(c1)],
            sum(c1[j][i] * x[i] for i in findnz(c1[j])[1]) <= 1 + sum(min.(findnz(c1[j])[2], 0))
        )
        #=
        for i in 1:length(c1)
            MOI.set(model_rebuild, Gurobi.ConstraintAttribute("Lazy"), con1[i], 3)
        end
       =#
    end
    
    if length(c2) > 0
        @constraint(
            model_rebuild,
            con2[j in 1:length(c2)],
            sum(c2[j][i] * x[i] for i in findnz(c2[j])[1]) <= 1 + sum(min.(findnz(c2[j])[2], 0))
        )
        if l3 >= 8000
            for i in 1:length(c2)
                MOI.set(model_rebuild, Gurobi.ConstraintAttribute("Lazy"), con2[i], -1)
            end
        end
    end
    
    if length(c3) > 0
        @constraint(
            model_rebuild,
            con3[j in 1:length(c3)],
            sum(c3[j][i] * x[i] for i in findnz(c3[j])[1]) <= 1 + sum(min.(findnz(c3[j])[2], 0))
        )
        if l2 > 8000
            for i in 1:length(c3)
                MOI.set(model_rebuild, Gurobi.ConstraintAttribute("Lazy"), con3[i], -1)
            end
        end
    end
    
    if length(c4) > 0
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

function rebuild_model_org(con_matrix, con_ub, con_lb, var_ub, var_lb, var_type, obj_coef, sp, bin_to_org, is_min)
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
    model_rebuild = Model()#Gurobi.Optimizer)

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
    
    set_pack = copy(sp)
    set_pack = remap_cons(set_pack, bin_to_org, var_lb);
    @constraint(model_rebuild, 
        c0[j in 1:length(set_pack)],
        sum(set_pack[j][i] * x[i] for i in findnz(set_pack[j])[1]) <= 1
    )
    
    #set objective function
    (is_min) ? (@objective(model_rebuild, Min, sum(obj_coef[i] * x[i] for i in 1:var_num))) : (@objective(model_rebuild, Max, sum(obj_coef[i] * x[i] for i in 1:var_num)))
    
    
    return model_rebuild
end