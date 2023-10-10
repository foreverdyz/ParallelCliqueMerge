#constraint_process.jl

include("binary_list.jl")
include("convert_to_standard.jl")
include("standard_con_process.jl")


function process_all_cons(A::AbstractSparseArray, con_ub::AbstractVector, con_lb::AbstractVector, var_ub::AbstractVector, var_lb::AbstractVector, var_type::AbstractVector)
    #build dicts between original variables and binary variables
    org_to_bin, _ = binary_list(var_type)
    #get m constraints
    m = size(A)[1]

    binary_number = 0
    for i in var_type
        (i > 1) && (binary_number += 1)
    end

    #record all set_pack constraints, like x1 + x2 + x3 <= 1, then terms like [1,2,3]
    set_pack = Vector{Int64}[]
    #for knapsack constraint, (a::SparseVector, b::Real), where b is RHS
    knapsack_set = Tuple{Any, Any}[]
    #for constraints like x1 = 1., fix_set[1] = 1; or if x1 = 0, fix_set[1] = -1
    fix_set = spzeros(binary_number)
    #indicator for whether the problem is feasible
    is_feasible = true

    #use for-loop to operate all constraints
    for j in 1:m
        #first we convert a general constraint to a standard form, which only has binary terms
        #id implies number of standard constraints we get from the general constraint
        id, a, b = convert_to_standard(A[j,:], con_ub[j], con_lb[j], var_ub, var_lb, var_type, org_to_bin)
        if id == 1
            dropzeros!(a)
            #newid could 0 (problem infeasible), -1 (fix variables), 1 (knapsack con), and 2 (set pack con)
            newid, con, new_b = standard_con_process(a, b)
            if newid == 0
                println(j)
                is_feasible =  false
                break
            elseif newid == -1
                #for fixed variables
                for i in findnz(a)[1]
                    (a[i] > 0) ? (fix_set[i] = -1) : (fix_set[i] = 1)
                end
            elseif newid == 1
                push!(knapsack_set, (con, new_b))
            else
                (length(findnz(con)[1]) > 1) && (push!(set_pack, findnz(con)[1]))
            end
        elseif id == 2
            dropzeros!(a[1])
            dropzeros!(a[2])
            newid, con, new_b = standard_con_process(a[1], b[1])
            if newid == 0
                println(j)
                is_feasible =  false
                break
            elseif newid == -1
                for i in findnz(a[1])[1]
                    (a[1][i] > 0) ? (fix_set[i] = -1) : (fix_set[i] = 1)
                end
            elseif newid == 1
                push!(knapsack_set, (con, new_b))
            else
                push!(set_pack, findnz(con)[1])
            end

            newid, con, new_b = standard_con_process(a[2], b[2])
            if newid == 0
                is_feasible =  false
                break
            elseif newid == -1
                for i in findnz(a[2])[1]
                    (a[2][i] > 0) ? (fix_set[i] = -1) : (fix_set[i] = 1)
                end
            elseif newid == 1
                push!(knapsack_set, (con, new_b))
            else
                push!(set_pack, findnz(con)[1])
            end
        else
            continue
        end
    end
    #output: feasibility?, set_pack list, knapsack list, fix_set as a sparse vector
    (is_feasible) ? (return true, set_pack, knapsack_set, fix_set) : (return false, false, false, false)
end