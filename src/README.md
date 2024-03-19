# Description for all Scripts

## Folder info_prepare

This folder includes scripts extracting model information, e.g. constraints coefficients, variables' types and bounds, for sequel clique merging.

## Folder clq_merge

This folder includes scripts for clique detecting, conflict graph building, clique extending, and domination checking between cliques.

## package_check.jl

Yon can check whether you have installed all required packages by running this script:
```julia
julia> include("package_check.jl")
```

## presolve_runtimes.jl

You can get and compare runtimes between the serial method (1 thread) and the parallel method (n threads) by this scrip:
```
$ julia --threads 16   #start julia with 16 threads (n = 16)
```
```julia
julia> include("presolve_runtimes.jl")
[ Info: You are using 16 threads.
[ Info: Using " julia --threads number_threads " to change the number of threads.
julia> test_runtime("30n20b8.mps", nthreads(), 0.00001, 3000, 1_000_000, 100_000)   #"30n20b8.mps" is the file name of the model.
```


## rebuild_model.jl

This script will build the reduced model by JuMP based on information provided by the presolving method.

## model_runtimes.jl

This script will read the mode from the .mps file (by scripts from [info_prepare](/src/info_prepare)), conduct the clique merging method (by scripts from [clq_merge](/src/clq_merge)), rebuild the reduced model (by [rebuild_model.jl](/src/rebuild_model.jl)), and solve both the original model and the reduced model by Gurobi 10.0.1.
```julia
julia> include("model_runtimes.jl")
julia> model_runtimes("30n20b8.mps", nthreads(), 0.00001, 3000, 1_000_000, 100_000)   #"30n20b8.mps" is the file name of the model.
```

## org_runtimes.jl

This script uses Gurobi to solve the model from .mps file
```julia
julia> include("org_runtimes.jl")
julia> org_runtimes("30n20b8.mps")   #"30n20b8.mps" is the file name of the model.
```

## pre_compile_model.mps
pre_compile_model.mps is a toy model for Julia to compile all functions we need.
