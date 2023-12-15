#model_info.jl
using JuMP
using SparseArrays

function model_info(filename::String)
    #read_from_file() is from JuMP
    model = read_from_file(filename)

    #info for variables and constraints
    #for variabels, collect:
    #variable names
    var_name = all_variables(model)
    #number of variables
    var_num = length(var_name)
    #list to denote variables' types: 2 imples binary, 1 implies integer, and 0 implies continuous
    var_type = spzeros(var_num)
    #record all integeral variables (note that, somtimes both binary and integer variables are denoted as integeral variables)
    for index in 1:var_num
        (is_integer(variable_by_name(model, string(var_name[index])))) && (var_type[index] = 1)
        (is_binary(variable_by_name(model, string(var_name[index])))) && (var_type[index] = 2)
    end

    #number of constraints, "include_variable_in_set_constraints=false" remove boxed constraints
    con_num = length(all_constraints(model, include_variable_in_set_constraints=false))
    #list to denote constraints' types: 1 imples >=, -1 implies <=, 0 implies ==
    con_type = zeros(con_num)

    #to get obj coefficient
    obj = MOI.get(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
    obj_coef = spzeros(var_num)
    for term in obj.terms
        obj_coef[term.variable.value] = term.coefficient
    end

    #to collect these necessary info, we need to relax the ip to an lp 
    undo = relax_integrality(model)
    #use JuMP._standard_form_matrix() to get info, this function is from JuMP
    #var_dic is a dictionary from varaible name => variable order
    #lb_list includes all lower bounds for variables + constraints
    #ub_list includes all upper bounds for varaibles + constraints
    #full_A is a matrix for all constraints and varaibles + slack variables
    _, lb_list, ub_list, full_A, _, _ = JuMP._standard_form_matrix(model)

    #extract info we need
    #split lb_list, ub_list into variables' and constraints'
    var_lb, con_lb = lb_list[1:var_num], lb_list[var_num+1:end]
    var_ub, con_ub = ub_list[1:var_num], ub_list[var_num+1:end]
    #get constraints' matrix only for original variables (excluding slack variables)
    con_matrix = full_A[1:con_num, 1:var_num]

    #denote binary variable:
    #if x is integer (var_type is 1) and 0<=x<=1, x is binary (change var_type to 2)
    for index in 1:var_num
        if var_type[index] == 1
            if var_lb[index] == 0 && var_ub[index] == 1
                var_type[index] = 2
            end
        end
    end

    #record objective sense
    if objective_sense(model) == MIN_SENSE
        is_min = true
    else
        is_min = false
    end

    #output info we need
    return con_matrix, con_ub, con_lb, var_ub, var_lb, var_type, obj_coef, is_min
end