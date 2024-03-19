#org_runtimes.jl

ENV["GRB_LICENSE_FILE"] = "/usr/local/gurobi/10.0.1/gurobi.lic"
using Gurobi
using BangBang
using JuMP

function org_runtimes(filename::String)
    model = read_from_file(filename)
    set_optimizer(model, Gurobi.Optimizer);
    #set_silent(model)
    set_optimizer_attribute(model, "LogFile", "res/logfile_org/"*filename[9:end-4]*"log.txt")
    set_attribute(model, "Threads", 10);
    set_time_limit_sec(model, 3600.0);
    optimize!(model)
    if has_values(model)
        return filename[9:end-4], objective_value(model), solve_time(model), node_count(model)
    else
        return filename[9:end-4], "Inf", solve_time(model), node_count(model)
    end
end