# Description for all Scripts

## Folder info_prepare

This folder includes scripts extracting model information, e.g. constraints coefficients, variables' types and bounds, for sequel clique merging.

## Folder presolve_cliuqes

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
julia> clq_merge_runtime("30n20b8.mps")   #"30n20b8.mps" is the file name of the model.
```

## presolve_runtimes_heuristic.jl

Similar to presolve_runtimes.jl

## rebuild_model.jl

This script will read the mode from the .mps file (by scripts from [info_prepare](/src/info_prepare)), conduct the clique merging method (by scripts from [presolve_cliques](/src/presolve_cliques)), and rebuild the reduced model.
