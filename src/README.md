The code includes three main functions

# Presolving:
Using the following code to presolve instances with our method, where the format of the filename "xxx.mps" r "xxx.mps.gz".
```julia
> include("presolve.jl")
> filename, presolve_time, _ = presolve(filename, threadnum) #threadnum is the number of threads you want to use
```
and the reduced model will be saved to the presolved_instances folder.

# Gurobi Solver Time Comparison
Use the following code to compare the Gurobi solver runtimes between the original and reduced models. threadnum is the number of threads to resolve,
random_seed is the random seed for the solver, and solver_thread is the number of threads for Gurobi.
```julia
> include(".jl")
> filename, org_time, org_node, reduced_time, reduced_node, presolve_time = test_runtime(filename, threadnum, random_seed, solver_thread)
```
and the logfile will be saved to the logfile folder.

# Other Solver
Since Julia JuMP cannot write the mps file with lazy constraints. You can use the following code to generate a reduced instance first
```julia
> include("presolve_other_solvers.jl")
> filename, presolve_time, _ = presolve_other_solvers(filename, threadnum) #threadnum is the number of threads you want to use
```
and the reduced model will be saved to presolved_instances_other_solvers. In the reduced model, all lazy cuts are constraints whose name
begins with "lazy_". You can select all of them and set them as lazy cut or user cut.
