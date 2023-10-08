#con_check.jl

function _is_feasible_or_set_pack(a::SparseVector, b::Real)
    s = 0
    (b == 1) ? (is_set_pack = true) : (is_set_pack = false)
    for i in findnz(a)[2]
        s += min(0, i)
        if is_set_pack
            (i != 1) && (is_set_pack = false)
        end
    end
    (s < b + 0.00001) ? (return true, is_set_pack) : (return false, false)
end
