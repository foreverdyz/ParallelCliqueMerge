#!bin/bash

julia --threads 1 presolve.jl 1
julia --threads 2 presolve.jl 2
julia --threads 4 presolve.jl 4
julia --threads 8 presolve.jl 8
julia --threads 16 presolve.jl 16
julia --threads 12 presolve.jl 12
