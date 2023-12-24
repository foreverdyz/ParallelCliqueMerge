# ParallelCliqueMerge

This archive is under the [MIT License](LICENSE).

The software and data in this repository are a snapshot of the software and data that were used in the research reported on in [Parallelized Conflict Graph Cut Generation](https://arxiv.org/abs/2311.03706) by Yongzheng Dai and Chen Chen.

## Language

The code is written in the [Julia programming language](https://julialang.org). Please visit the website for an installation guide. 

## Required packages 

To run the content of scripts in [src](/src), the following Julia packages are required: Gurobi, ThreadsX, BangBang, Random, SparseArrays, JuMP.

Yon can check whether you have installed these required packages by running [package_check.jl (in src)](/src/packge_check.jl)
```julia
julia> include("package_check.jl")
```

To install a package simply run

```julia
pkg> add PACKAGE_NAME    # Press ']' to enter the Pkg REPL mode.
```

## Tutorial

See [README.md (in src)](/src/README.md) for tutorial and examples.
