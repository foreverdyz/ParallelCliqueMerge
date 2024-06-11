# Description for all Scripts

## Language and Packages

Code is implemented in Julia. Required packages: "JuMP", "Gurobi", "BangBang", "SparseArrays", "Random"

We have added these packages to MPC server. The only thing we want to mention here is that:

we are using an alpha version of JuMP, and the prototype can be got by:

```julia
pkg> add JuMP  #Press ']' to enter the Pkg REPL mode.

pkg> add MathOptInterface#master
``` 


## Source Code

### Folder: mipdata

This folder includes all mps file from MIPLIB 2017 Benchmark Collection.

### Folder: res

It is an empty folder to cache exp results.

### Folder: info_prepare

This folder includes scripts to read mps file, one-round domain propagation, detect binary variables (if needed), cliques, and knapsacks.

Generally, using

```julia
julia> include("info_prepare/generate_info.jl");

julia> info = generate_info(filename);
```

can get necessary info for our method.

### Folder: clq_merge

This folder includes scripts to generate conflict graph cut in parallel (or serial): detect cliques from knapsacks, construct conflict graph, extend cliques, and merge cliques. 

Generally, using: 

```julia
julia> include("clq_merge/generate_sets.jl")
```

to get function "generate_sets()", which can return conflict graph cuts based on info from "generate_info()".

### model_build.jl

This script build MIP model by JuMP based on the original model and our conflict graph cuts. Note that, this model cannot be written to file since it includes lazy cuts (featured by Gurobi).

### get_reduced_model.jl

This script returns the reduced model for a provided model name.

### Other Files

All other files are built for our experiments, which are introduced in later session.

## Experiments

### Parallel Performance

First, start Julia with number of threads you want to test the parallel performance, e.g. 8 threads:

```
$ julia --threads 8
```

Then get parallel cut generation time:

```julia
julia> include("presolve_runtime.jl")
```

Note that, this scripts will run all cases, which is time-consuming. We suggest to run with following two ways:

```
$ nohup Julia --threads 8 presolve_runtime.jl > res/presolve_runtimes_8_threads.txt
```

or change our code a little bit:

```
$ vi presolve_runtime.jl
``` 

in line 57, for-loop, change "list" to "list[1:10]" or whatever cases you want to test.

Additionally, we use macro "@time", which includes cg time and compile time. You may need to remove them manually.

### Solver Time Test

First, you may also want to start Julia with multiple threads (but unnecessary):

```
$ julia --threads 16
```

Then get solver time and nodes number (both org model and our reduced model) and collect them in "res/exp_cmp_res.csv" by:

```julia
julia> include("test_cmp_runtimes.jl")
```

In the csv file, four columns are "solver_time_org, node_count_org, solver_time_reduced, node_count_reduced".


Note that, this scripts will run all cases, which is very time-consuming. We suggest to run with following two ways:

```
$ nohup julia --threads 16 presolve_runtime.jl
```

or change our code a little bit:

```
$ vi test_cmp_runtimes.jl
``` 

in line 12, for-loop, change "list" to "list[1:10]" or whatever cases you want to test.
