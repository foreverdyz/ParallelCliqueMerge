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
[ Info: You are using 6 threads.
[ Info: Using " julia --threads number_threads " to change the number of threads.
```
