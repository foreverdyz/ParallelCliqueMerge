#!/bin/bash

for ((i=1; i<=99; i++))
do
    julia --threads 12 test_runtime.jl 12 0 $i $i
done

