#parameters.jl
using JuMP, Gurobi
using NLPModelsJuMP
using SparseArrays

model = read_from_file("parallel_mip/miplib/graphdraw-domain.mps")

var_name = all_variables(model)

var_num = length(var_name)
con_num = length(all_constraints(model, include_variable_in_set_constraints=false))

var_type = zeros(var_num)

for index in 1:var_num
    (is_integer(variable_by_name(model, string(var_name[index])))) && (var_type[index] = 1)
end

undo = relax_integrality(model);

var_dic, lb_list, ub_list, full_A, _, _ = JuMP._standard_form_matrix(model)

var_lb, con_lb = lb_list[1:var_num], lb_list[var_num+1 : end]
var_ub, con_ub = ub_list[1:var_num], ub_list[var_num+1:end]
con_matrix = full_A[1:con_num, 1:var_num]

for index in 1:var_num
    if var_type[index] > 0
        if var_lb[index] == 0 && var_ub[index] == 1
            var_type[index] = 2
        end
    end
end

con_type = zeros(con_num)

for index in 1:con_num
    if con_lb[index] == -Inf
        con_type[index] = -1
    elseif con_ub[index] == Inf
        con_type[index] = 1
    end
end

model_rebuild = Model(Gurobi.Optimizer)
@variable(model_rebuild, var_ub[i] >= x[i in 1:var_num] >= var_lb[i])
for i in 1:var_num
    (var_type[i] == 2) && (set_binary(x[i]))
    (var_type[i] == 1) && (set_integer(x[i]))
end

for j in 1:con_num
    if con_type[j] == 1
        @constraint(model_rebuild, sum(con_matrix[j, i] * x[i] for i in findnz(con_matrix[j, :])[1]) >= con_lb[j])
    elseif con_type[j] == -1
        @constraint(model_rebuild, sum(con_matrix[j, i] * x[i] for i in findnz(con_matrix[j, :])[1]) <= con_ub[j])
    else
        @constraint(model_rebuild, sum(con_matrix[j, i] * x[i] for i in findnz(con_matrix[j, :])[1]) == con_lb[j])
    end
end


nlp = MathOptNLPModel(model)
tmp = zeros(nlp.meta.nvar)
obj = NLPModelsJuMP.grad(nlp, tmp)
if objective_sense(model) == MIN_SENSE
    @objective(model_rebuild, Min, sum(obj[i] * x[i] for i in 1:var_num))
else
    @objective(model_rebuild, Max, sum(obj[i] * x[i] for i in 1:var_num))
end

undo;

set_optimizer(model, Gurobi.Optimizer);