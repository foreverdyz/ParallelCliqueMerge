#model_info.jl

using JuMP #use it to read the model
using SparseArrays #constraints info are cached in sparse form

"""
    model_info(filename::String)
This function read the model from the file (named filename), and extract some necessary info for our sequel presolve
"""
function model_info(filename::String)
    model = read_from_file(filename)

    #info for variables and constraints
    #for variabels, collect:
    #variable names
    var_name = all_variables(model)
    #number of variables
    var_num = length(var_name)
    #list to denote variables' types: 2 imples binary, 1 implies integer, and 0 implies continuous
    var_type = spzeros(var_num)
    #record all integeral variables (note that, both binary and integer variables are denoted as integeral variables)
    for index in 1:var_num
        (is_integer(variable_by_name(model, string(var_name[index])))) && (var_type[index] = 1)
        (is_binary(variable_by_name(model, string(var_name[index])))) && (var_type[index] = 2)
    end

    #number of constraints
    con_num = length(all_constraints(model, include_variable_in_set_constraints=false))
    #list to denote constraints' types: 1 imples >=, -1 implies <=, 0 implies ==
    con_type = zeros(con_num)
    
    #to get obj coefficient
    obj = MOI.get(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}());
    obj_coef = spzeros(var_num);
    for term in obj.terms
        obj_coef[term.variable.value] = term.coefficient
    end
    
    #to collect these necessary info, we need to relax the ip to an lp 
    undo = relax_integrality(model)
    #use JuMP._standard_form_matrix to get info
    #var_dic is a dictionary from varaible name => variable order
    #lb_list includes all lower bounds for variables + contraints
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
    #if x is integer and 0<=x<=1, x is binary
    for index in 1:var_num
        if var_type[index] > 0
            if var_lb[index] == 0 && var_ub[index] == 1
                var_type[index] = 2
            end
        end
    end

    #determine whether it is a min or max problem
    if objective_sense(model) == MIN_SENSE
        is_min = true
    else
        is_min = false
    end

    #here we use another way to store matrix, since using sparse matrix is slow in following code 
    con_set = Vector{Int64}[]
    con_coef = []
    for i in 1:size(con_matrix)[1]
        local ind, coef = findnz(con_matrix[i, :])
        push!(con_set, ind)
        push!(con_coef, coef)
    end
    
    return con_matrix, con_set, con_coef, con_ub, con_lb, var_ub, var_lb, var_type, obj_coef, is_min
end