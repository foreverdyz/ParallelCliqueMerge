The code includes three main functions

# Presolving:
```julia
> include("presolve.jl")
> presolve(filename)
``` 

# Gurobi Solving
Run:
~$free -g
to check whether have enough memory again. Generally our results should be reproduced with more than 14GB in Mem.
Run:
~$nohup bash randomseed\_0.sh
to use gurobi solving both original and reduced models. Results are saved to the folder res.
We want to emphasize several points:
1, we can complete all experiments by run randomseed\_0.sh within one round in most situations. If sometimes the program is killed (generally due to memory issue because this server shares memory with other accounts), you can check res/runtime\_1\_9912.csv and found the last line number (e.g. line 66). Then open randomseed\_0.sh and change the start case in the for-loop from 1 to 67 and rerun the randomseed\_0.sh.
2, if you want to change the number of threads, just change julia --threads 12 to whatever number of threads yuo want. But threads 16 will make the experiments broken in some cases, including original and reduced instances.
3, if you want to change the random seed, just change julia --threads 12 test\_runtime.jl 12 0 from 0 to whatever seed you want.
4, We run the experiments multiple times, we found in most runs, the runtime comparison is 91.5%; some runs are around 95% or 85% but we do not know why. If you cannot reporduce the 91.5% res, you may 1, rerun both presolve.sh and randomseed\_0.sh in a different time; 2, you can also try our saved res: cp presolved\_data\_12\_author\_run to presolved\_data\_12 and cp presolve\_runtime\_12\_author\_run.csv to presolve\_runtime\_12.csv and run randomseed\_0.sh again.





