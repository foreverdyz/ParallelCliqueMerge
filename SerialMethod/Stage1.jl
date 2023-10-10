#Stage1.jl

#load model information 
include("model_info.jl")

#presolve all variables
include("variable_bound.jl")
con_matrix, con_ub, con_lb, var_ub, var_lb, var_type = singleton_update(con_matrix, con_ub, con_lb, var_ub, var_lb, var_type);
var_ub, var_lb = bound_strengthen(con_matrix, con_ub, con_lb, var_ub, var_lb, var_type);

#presolve all constraints
include("constraint_process.jl")
fes, set_pack, knapsack_set, fix_set = process_all_cons(con_matrix, con_ub, con_lb, var_ub, var_lb, var_type);

if fes == false
    println(fes)
else
    if length(findnz(fix_set)[1]) == 0
        println("no fixed variable")
    else
        println("There are ", length(findnz(fix_set)[1]), " fixed varaibles")
    end
    if length(set_pack) == 0
        println("no set pack constraints")
    else
        println("There are ", length(set_pack), " set pack constraints")
    end
    if length(knapsack_set) == 0
        println("no knapsack constriant")
    else
        I = []
        for i in 1:length(knapsack_set)
            local x = findnz(knapsack_set[i][1])[2]
            if length(x) > 1
                x_sort = sort(x, rev = true)
                (x_sort[1] + x_sort[2] > knapsack_set[i][2]) && (push!(I, i))
            end
        end
        if I == []
            println("No useful knapsack set")
        else
            println("There are ", length(I), " useful knapsack set")
        end
    end
end