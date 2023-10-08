#con_check.jl

function _is_feasible(a::SparseVector, b::Real)
    s = 0
    for i in findnz(a)[2]
        s += min(0, i)
    end
    (s < b + 0.00001) ? (return true) : (return false)
end