#standard_con_process.jl

function standard_con_process(a::SparseVector, b::Real)
    #check feasibility
    id = _is_feasible(a, b)
    if id == 1
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


function _is_feasible(a::SparseVector, b::Real)
    #record the min of LHS
    s = 0
    for i in findnz(a)[2]
        s += min(0, i)
    end
    if s > b
        #0 implies this problem is infeasible
        return 0
    elseif s < b + 0.0001 && s > b - 0.0001
        #-1 implies all variables in this constraint should be fixed to the lower bound (depends on a[i])
        return -1
    else
        #1 implies a knapsack constraint
        return 1
    end
end

function _expand_constraint(a::SparseVector, b::Real)
    binary_length = length(a)
    #build an empty expended constraint
    con = spzeros(2*binary_length)
    is_set_pack = true
    for i in findnz(a)[1]
        local c = a[i]
        if c > 0
            con[i] = c
        else
            con[i + binary_length] = -c
            b -= c
        end
        if is_set_pack
            (abs(c) != 1) && (is_set_pack = false)
        end
    end
    if is_set_pack
        (b != 1) && (is_set_pack = false)
    end
    return con, b, is_set_pack
end