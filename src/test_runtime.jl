#diff_runtime.jl

using Gurobi
using JuMP
using DelimitedFiles

include("presolve.jl")

path3 = "presolved_instances/";
mkpath("logfile")


function(filename, threadnum, random_seed, solver_thread)
    
    #for orginal model
    model_0 = read_from_file(filename);
    set_optimizer(model_0, Gurobi.Optimizer);
    set_time_limit_sec(model_0, 3600.0);
    set_optimizer_attribute(model_0, "LogFile", "logfile/"*filename[1:end-4]*"_0log.txt");
    set_optimizer_attribute(model_0, "Threads", solver_thread);
    set_optimizer_attribute(model_0, "Seed", random_seed);
    optimize!(model_0);
    t0, n0 = solve_time(model_0), node_count(model_0);
    GC.gc();
    
    _, pre_time, tol, _, c1, c2, c3, c4, _, _, _, _ = presolve(filename, threadnum)
    #judge user cuts
    alpha = 1
    tol, c2, c3, c4 = presolve_info[iter, [3, 6, 7, 8]];
    iden = 0;
    if (c2+c3+c4) <= alpha*tol
	    iden = 1234;
    else
	    if (c2+c3) < alpha*tol
	        iden = 123;
	    elseif c2 < alpha*tol
	        iden = 12;
	    elseif c3 < alpha*tol
	        iden = 13;
    	else
	        iden = 1
	    end
    end
    
    if iden == 1
    #add nothing, all new cuts to user cuts pool but just original setpacking constriants (extended)
    model_1 = read_from_file(path3*filename);
    con = all_constraints(model_1; include_variable_in_set_constraints = false);
    for x in con
        if name(x)[1:2] == "la"
            #lazy_c1 -> root node, and all others in user cuts pool
            if string(name(x)[7]) != "1" || string( name(x)[8]) == "2" 
                MOI.set(model_1, Gurobi.ConstraintAttribute("Lazy"), x, -1)
            end
        end
    end
    set_optimizer(model_1, Gurobi.Optimizer);
    set_time_limit_sec(model_1, 3600.0);
    set_optimizer_attribute(model_1, "LogFile", "logfile/"*filename[1:end-4]*"_1log.txt");
    set_optimizer_attribute(model_1, "Seed", random_seed);
    set_optimizer_attribute(model_1, "Threads", solver_thread);
    optimize!(model_1);
    t1, n1 = solve_time(model_1), node_count(model_1);
    #GC.gc();
    end
    
    if iden == 12
    #add c1, c2
    model_2 = read_from_file(path3*filename);
    con = all_constraints(model_2; include_variable_in_set_constraints = false);
    count_temp_2 = 0;
    for x in con
        if name(x)[1:2] == "la"
            #lazy_c1 and lazy_c2 -> root node, and all others in user cuts pool
            if string(name(x)[7] )!= "1" && string(name(x)[7]) != "2"  
                MOI.set(model_2, Gurobi.ConstraintAttribute("Lazy"), x, -1)
            elseif string(name(x)[8]) == "2"
                MOI.set(model_2, Gurobi.ConstraintAttribute("Lazy"), x, -1)
            else
                if string(name(x)[7]) == "2"  
                    count_temp_2 += 1;
                end
            end
        end
    end
    set_optimizer(model_2, Gurobi.Optimizer);
    set_time_limit_sec(model_2, 3600.0);
    set_optimizer_attribute(model_2, "LogFile", "logfile/"*filename[1:end-4]*"_2log.txt");
    set_optimizer_attribute(model_2, "Threads", solver_thread);
    set_optimizer_attribute(model_2, "Seed", random_seed);
    optimize!(model_2);
    t1, n1 = solve_time(model_2), node_count(model_2);
    #GC.gc();
    end
   
    if iden == 13
    #add c1, c3
    model_3 = read_from_file(path3*filename);
    con = all_constraints(model_3; include_variable_in_set_constraints = false);
    count_temp_3 = 0;
    for x in con
        if name(x)[1:2] == "la"
            #lazy_c1 and lazy_c3 -> root node, and all others in user cuts pool
            if string(name(x)[7]) != "1" && string(name(x)[7]) != "3" 
                MOI.set(model_3, Gurobi.ConstraintAttribute("Lazy"), x, -1)
            elseif string(name(x)[8]) == "2" 
                MOI.set(model_3, Gurobi.ConstraintAttribute("Lazy"), x, -1)
            else
                if string(name(x)[7]) == "3"  
                    count_temp_3 += 1;
                end
            end
        end
    end
    set_optimizer(model_3, Gurobi.Optimizer);
    set_time_limit_sec(model_3, 3600.0);
    set_optimizer_attribute(model_3, "LogFile", "logfile/"*filename[1:end-4]*"_3log.txt");
    set_optimizer_attribute(model_3, "Threads", solver_thread);
    set_optimizer_attribute(model_3, "Seed", random_seed);
    optimize!(model_3);
    t1, n1 = solve_time(model_3), node_count(model_3);
    #GC.gc();
    end

    if iden == 123
    #add c1, c2, c3
    model_4 = read_from_file(path3*filename);
    con = all_constraints(model_4; include_variable_in_set_constraints = false);
    for x in con
        if name(x)[1:2] == "la"
            #lazy_c1, lazy_c2, and lazy_c3 -> root node, and all others in user cuts pool
            if string(name(x)[7]) == "4"
                MOI.set(model_4, Gurobi.ConstraintAttribute("Lazy"), x, -1)
            elseif string(name(x)[8]) == "2" 
                MOI.set(model_4, Gurobi.ConstraintAttribute("Lazy"), x, -1)
            end
        end
    end
    set_optimizer(model_4, Gurobi.Optimizer);
    set_time_limit_sec(model_4, 3600.0);
    set_optimizer_attribute(model_4, "LogFile", "logfile/"*filename[1:end-4]*"_4log.txt");
    set_optimizer_attribute(model_4, "Threads", solver_thread);
    set_optimizer_attribute(model_4, "Seed", random_seed);
    optimize!(model_4);
    t1, n1 = solve_time(model_4), node_count(model_4);
    #GC.gc();
    end
    
    if iden == 1234
    #add c1, c2, c3, c4
    model_5 = read_from_file(path3*filename);
    con = all_constraints(model_5; include_variable_in_set_constraints = false);
    count_temp_5 = 0;
    for x in con
        if name(x)[1:2] == "la"
            #lazy_c1, lazy_c2,  lazy_c3, and lazy_c4 -> root node, and all others in user cuts pool
            if string(name(x)[8]) == "2" 
                MOI.set(model_5, Gurobi.ConstraintAttribute("Lazy"), x, -1)
            end
            if string(name(x)[7]) == "4"
                count_temp_5 += 1;
            end
        end
    end
    set_optimizer(model_5, Gurobi.Optimizer);
    set_time_limit_sec(model_5, 3600.0);
    set_optimizer_attribute(model_5, "LogFile", "logfile/"*filename[1:end-4]*"_5log.txt");
    set_optimizer_attribute(model_5, "Threads", solver_thread);
    set_optimizer_attribute(model_5, "Seed", random_seed);
    optimize!(model_5);
    t1, n1 = solve_time(model_5), node_count(model_5);
    #GC.gc();
    end  
    return filename, t0, n0, t1, n1, pre_time
end
