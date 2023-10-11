using JuMP
#using Gurobi
using SparseArrays

function model_info(filename::String)
    model = read_from_file("parallel_mip/miplib/mipdata1-25/"*filename*".mps")

    #info for variables and constraints
    #for variabels, collect:
    #variable names
    var_name = all_variables(model)
    #number of variables
    var_num = length(var_name)
    #list to denote variables' types: 2 imples binary, 1 implies integer, and 0 implies continuous
    var_type = zeros(var_num)
    #record all integeral variables (note that, both binary and integer variables are denoted as integeral variables)
    for index in 1:var_num
        (is_integer(variable_by_name(model, string(var_name[index])))) && (var_type[index] = 1)
        (is_binary(variable_by_name(model, string(var_name[index])))) && (var_type[index] = 2)
    end

    #number of constraints
    con_num = length(all_constraints(model, include_variable_in_set_constraints=false))
    #list to denote constraints' types: 1 imples >=, -1 implies <=, 0 implies ==
    con_type = zeros(con_num)

    #info of objective function
    obj_coef = zeros(var_num)
    obj_abs = MOI.get(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
    for (i, term) in enumerate(obj_abs.terms)
        obj_coef[i] = term.coefficient
    end

    #to collect these necessary info, we need to relax the ip to an lp 
    undo = relax_integrality(model)

    #use JuMP._standard_form_matrix to get info
    #var_dic is a dictionary from varaible name => variable order
    #lb_list includes all lower bounds for variables + contraints
    #ub_list includes all upper bounds for varaibles + constraints
    #full_A is a matrix for all constraints and varaibles + slack variables
    var_dic, lb_list, ub_list, full_A, _, _ = JuMP._standard_form_matrix(model)

    #extract info we need
    #split lb_list, ub_list into variables' and constraints'
    var_lb, con_lb = lb_list[1:var_num], lb_list[var_num+1:end]
    var_ub, con_ub = ub_list[1:var_num], ub_list[var_num+1:end]
    #get constraints' matrix only for original variables (excluding slack variables)
    con_matrix = full_A[1:con_num, 1:var_num]

    #denote binary variable:
    #if x is integer and 0<=x<=1, x is binary
    for index in 1:var_num
        if var_type[index] > 0
            if var_lb[index] >= 0 && var_ub[index] <= 1
                var_type[index] = 2
            end
        end
    end
    return con_matrix, con_ub, con_lb, var_ub, var_lb, var_type
end

#=
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
#now we rebuild the model:
#check we get correct info of model
model_rebuild = Model(Gurobi.Optimizer)

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
    else
        @constraint(model_rebuild, sum(con_matrix[j, i] * x[i] for i in findnz(con_matrix[j, :])[1]) == con_lb[j])
    end
end

#set objective function
if objective_sense(model) == MIN_SENSE
    @objective(model_rebuild, Min, sum(obj_coef[i] * x[i] for i in 1:var_num))
else
    @objective(model_rebuild, Max, sum(obj_coef[i] * x[i] for i in 1:var_num))
end

#relax model for check our rebuilding model
undo();

set_optimizer(model, Gurobi.Optimizer);
#=
optimize!(model_rebuild)
optimize!(model)

@show objective_value(model)
@show objective_value(model_rebuild)
=#
=#