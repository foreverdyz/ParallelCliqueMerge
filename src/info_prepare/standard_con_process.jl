#standard_con_process.jl

using SparseArrays

function standard_con_process(a::SparseVector, b::Real)
    #check feasibility, id could be 0 (infeasible), 1 (feasible), -1 (all variables need to be fixed lower bounds)
    id = _is_feasible(a, b)
    if id == 1
        #transfer the pure binary constraint to the form with complementary variables
        con, newb, newid = _expand_constraint(a, b)
        #1 imples knapsack constraint, and 2 implies set pack constraint
        (newid == false) ? (return 1, con, newb) : (return 2, con, newb)
    else
        #here id might be -1 or 0
        #0 implies this problem is infeasible
        #-1 imples variables should be fixed to the lower bound (depends on a[i])
        return id, false, false
    end
end

#this function tries to check whether this constraint is feasible (1), infeasible (0), and all variables need to fixed to the lower bound (-1)
function _is_feasible(a::SparseVector, b::Real)
    #record the min of LHS
    s = 0
    #ax <= b with only pure binary variables
    #if a[i] > 0, x = 0 to min
    #if a[i] < 0, x = 1 to min
    for i in findnz(a)[2]
        s += min(0, i)
    end
    #s > b => min(ax) > b
    if s > b
        #0 implies this problem is infeasible
        return 0
    #s = b => min(ax) = b
    elseif s < b + 0.0001 && s > b - 0.0001
        #-1 implies all variables in this constraint should be fixed to the lower bound (depends on a[i])
        return -1
    else
        #1 implies a knapsack constraint
        return 1
    end
end

#this function tries to add complementary variables to the pure binary constraint
#and also tries to check whether the constraint is a set-pack
function _expand_constraint(a::SparseVector, b::Real)
    binary_length = length(a)
    #build an empty expended constraint
    con = spzeros(2*binary_length)
    is_set_pack = true
    for i in findnz(a)[1]
        local c = a[i]
        #add complementary varaible when c < 0
        if c > 0
            con[i] = c
        else
            #-c > 0 
            con[i + binary_length] = -c
            #cx => c(1-x), and move c to RHS
            b -= c
        end
        if is_set_pack
            (abs(c) != 1) && (is_set_pack = false)
        end
    end
    if is_set_pack
        #RHS of set-pack is 1 
        (b != 1) && (is_set_pack = false)
    end
    #output
    return con, b, is_set_pack
end