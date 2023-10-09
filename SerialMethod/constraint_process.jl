#constraint_process.jl

include("binary_list.jl")
include("convert_to_standard.jl")
include("standard_con_process.jl")


function process_all_cons(A::AbstractSparseArray, con_ub::AbstractVector, con_lb::AbstractVector, var_ub::AbstractVector, var_lb::AbstractVector, var_type::AbstractVector)
    #build dicts between original variables and binary variables
    org_to_bin, bin_to_org = binary_list(var_type)
    #get m constraints and n variables
    m, n = size(A)

    set_pack = Vector{Int64}[]
    #for knapsack constraint, (a::SparseVector, b::Real)
    knapsack_set = Tuple{Any, Any}[]
    #for constraints like x1 = 1., element likes (1 ,1) where (1,.) for x1, and (.,1) for rhs 1
    fix_set = Tuple{Int64, Int64}[]
    is_feasible = true
    #use for-loop to operate all constraints
    for j in 1:m
        id, a, b = convert_to_standard(A[j,:], con_ub[j], con_lb[j], var_ub, var_lb, var_type, org_to_bin)
        if id == 1
            dropzeros!(a)
            newid, con, new_b = standard_con_process(a, b)
            if newid == 0
                is_feasible =  false
                break
            elseif newid == -1
                for i in findnz(a)[1]
                    (a[i] > 0) ? (push!(fix_set, (i, 0))) : (push!(fix_set, (i, 1)))
                end
            elseif newid == 1
                push!(knapsack_set, (con, new_b))
            else
                push!(set_pack, findnz(con)[1])
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
                for i in findnz(a)[1]
                    (a[i] > 0) ? (push!(fix_set, (i, 0))) : (push!(fix_set, (i, 1)))
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
                for i in findnz(a)[1]
                    (a[i] > 0) ? (push!(fix_set, (i, 0))) : (push!(fix_set, (i, 1)))
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
    (is_feasible) ? (return true, set_pack, knapsack_set, fix_set) : (return false, false, false, false)
end