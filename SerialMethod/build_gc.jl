#build_gc.jl

"""
Here we need number of binary variables
"""
function build_gc(set_pack::Vector{Vector{Int64}}, main_cliques::Vector{Vector{Int64}}, clique_set::Vector{Vector{Int64}}, binary_number::Int64)
    gc = _init_gc(binary_number)
end

function _init_gc(binary_number::Int64)
    gc = spzeros(Int64, 2*binary_number, 2*binary_number)
end