#constraint_process.jl

include("convert_to_standard.jl")
include("standard_con_process.jl")

#for all constraints, transform them to pure binary form with complementary variables
function process_all_cons(A::AbstractSparseArray, con_ub::AbstractVector, con_lb::AbstractVector, var_ub::AbstractVector, var_lb::AbstractVector, var_type::AbstractVector, org_to_bin::Dict{Int64, Int64}, bin_to_org::Dict{Int64, Int64})
    #get m constraints
    m = size(A)[1]

    binary_number = length(org_to_bin)

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
        #id imples this constraint provide 0, 1, or 2 standard constraints
        if id == 1
            dropzeros!(a)
            #newid could be 0 (problem infeasible), -1 (fix variables), 1 (knapsack con), and 2 (set pack con)
            newid, con, new_b = standard_con_process(a, b)
            if newid == 0
                println(j)
                #the problem is infeasible
                is_feasible =  false
                break
            elseif newid == -1
                #for fixed variables
                for i in findnz(a)[1]
                    (a[i] > 0) ? (fix_set[i] = -1) : (fix_set[i] = 1)
                end
            elseif newid == 1
                local x = findnz(con)[2]
                if length(x) > 1
                    #check whether the summation of the two largest coefficients is greater than RHS or not
                    x_sort = sort(x, rev=true)
                    (x_sort[1] + x_sort[2] > new_b) && (push!(knapsack_set, (con, new_b)))
                end
            else
                (length(findnz(con)[1]) > 1) && (push!(set_pack, findnz(con)[1]))
            end
        elseif id == 2
            dropzeros!(a[1])
            dropzeros!(a[2])
            #the first constraint
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
                local x = findnz(con)[2]
                if length(x) > 1
                    x_sort = sort(x, rev=true)
                    (x_sort[1] + x_sort[2] > new_b) && (push!(knapsack_set, (con, new_b)))
                end
            else
                (length(findnz(con)[1]) > 1) && (push!(set_pack, findnz(con)[1]))
            end
            #the second constraint
            newid, con, new_b = standard_con_process(a[2], b[2])
            if newid == 0
                is_feasible =  false
                break
            elseif newid == -1
                for i in findnz(a[2])[1]
                    (a[2][i] > 0) ? (fix_set[i] = -1) : (fix_set[i] = 1)
                end
            elseif newid == 1
                local x = findnz(con)[2]
                if length(x) > 1
                    x_sort = sort(x, rev=true)
                    (x_sort[1] + x_sort[2] > new_b) && (push!(knapsack_set, (con, new_b)))
                end
            else
                (length(findnz(con)[1]) > 1) && (push!(set_pack, findnz(con)[1]))
            end
        else
            continue
        end
    end
    for i in 1:binary_number
        j = bin_to_org[i]
        if var_lb[j] == var_ub[j]
            (var_lb[j] == 1) ? (fix_set[i] = 1) : (fix_set[i] = -1)
        end
    end
    #output: feasibility?, set_pack list, knapsack list, fix_set as a sparse vector
    (is_feasible) ? (return true, set_pack, knapsack_set, fix_set) : (return false, false, false, false)
end