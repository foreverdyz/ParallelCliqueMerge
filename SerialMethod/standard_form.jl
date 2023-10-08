#stand_form.jl

using SparseArrays

function convert_to_standard(a::SparseVector, ub::Real, lb::Real, var_ub::AbstractVector, var_lb::AbstractVector, var_type::AbstractVector, org_to_bin::Dict{Int64, Int64})
    if lb == -Inf
        con, b =  _con_ub_convert(a, ub, var_ub, var_lb, var_type, org_to_bin)
        (con == false) ? (return 0, false, false) : (return 1, con, b)
    elseif ub == Inf
        con, b =  _con_lb_convert(a, lb, var_ub, var_lb, var_type, org_to_bin)
        (con == false) ? (return 0, false, false) : (return 1, con, b)
    else
        (lb != ub) && (println("Warning: there is a constaint that lb and ub is not equal!!!!!!"))
        id, con, b =  _con_b_convert(a, ub, lb, var_ub, var_lb, var_type, org_to_bin)
        return id, con, b
    end
end

function _con_ub_convert(a::SparseVector, b::Real, var_ub::AbstractVector, var_lb::AbstractVector, var_type::AbstractVector, org_to_bin::Dict{Int64, Int64})
    #record number of binary variables in this constriant
    binary_number = 0
    #record lower bound of all non binary variables
    delta = 0
    #build a sparse vector only for binary variables
    con = spzeros(length(org_to_bin))

    #partially scan all non-zero coefficient variables
    for i in findnz(a)[1]
        if var_type[i] == 2
            binary_number += 1
            con[org_to_bin[i]] = a[i]
        else
            if a[i] >0
                #if this variable with lb inf, we cannot remove this variable from constraint; so cannot get a standard format
                (var_lb[i] == -Inf) ? (return false, false) : (delta += a[i]*var_lb[i])
            else
                (var_ub[i] == Inf) ? (return false, false) : (delta += a[i]*var_ub[i])
            end 
        end
    end

    #if no binary, do not have constraint we want
    (binary_number < 1) ? (return false, false) : (return con, b - delta)

end

function _con_lb_convert(a::SparseVector, b::Real, var_ub::AbstractVector, var_lb::AbstractVector, var_type::AbstractVector, org_to_bin::Dict{Int64, Int64})
    #record number of binary variables in this constriant
    binary_number = 0
    #record lower bound of all non binary variables
    delta = 0
    #build a sparse vector only for binary variables
    con = spzeros(length(org_to_bin))

    #partially scan all non-zero coefficient variables
    for i in findnz(a)[1]
        if var_type[i] == 2
            binary_number += 1
            con[org_to_bin[i]] = -a[i]
        else
            if a[i] >0
                #if this variable with lb inf, we cannot remove this variable from constraint; so cannot get a standard format
                (var_ub[i] == Inf) ? (return false, false) : (delta += a[i]*var_ub[i])
            else
                (var_lb[i] == -Inf) ? (return false, false) : (delta += a[i]*var_lb[i])
            end 
        end
    end
    
    #if no binary, do not have constraint we want
    (binary_number < 1) ? (return false, false) : (return con, -b + delta)
end

function _con_b_convert(a::SparseVector, ub::Real, lb::Real, var_ub::AbstractVector, var_lb::AbstractVector, var_type::AbstractVector, org_to_bin::Dict{Int64,Int64})
    con1, b1 = _con_ub_convert(a, ub, var_ub, var_lb, var_type, org_to_bin)
    con2, b2 = _con_lb_convert(a, lb, var_ub, var_lb, var_type, org_to_bin)
    if con1 != false && con2 != false
        return 2, [con1, con2], [b1, b2]
    elseif con1 != false && con2 == false
        return 1, con1, b1
    elseif con1 == false && con2 != false
        return 1, con2, b2
    else
        return 0, false, false
    end
end

