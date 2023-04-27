mutable struct TreeVectorAxes{N, NAME, LEAF}
    size::NTuple{N, Int}
    index::NTuple{N, UnitRange{Int}}
    leaf::LEAF
end
TreeVectorAxes() = TreeVectorAxes{0, (), Tuple{}}((), (), ())

function TreeVectorAxes(leaves, idx = 1)
    name = Vector{Symbol}()
    size = Vector{Int}()
    index = Vector{UnitRange{Int}}()
    leaf = Vector{Any}()
    for (n, l) in named_pairs(leaves)
        push!(name, n)
        if l isa Integer
            len = l
            push!(leaf, ())
        elseif l isa AbstractVector
            len = length(l)
            push!(leaf, ())
        else
            ll = TreeVectorAxes(l, idx)
            len = length(ll)
            push!(leaf, ll)
        end
        push!(size, len)
        push!(index, idx:(idx + len - 1))
        idx += len
    end
    _N = length(name)
    _NAME = (name...,)
    _size = (size...,)
    _index = (index...,)
    _leaf = (leaf...,)
    _LEAF = typeof(_leaf)
    TreeVectorAxes{_N, _NAME, _LEAF}(_size, _index, _leaf)
end

Base.isempty(::TreeVectorAxes) = false
Base.isempty(::TreeVectorAxes{0}) = true
Base.length(x::TreeVectorAxes) = sum(x.size)

named_pairs(x::NamedTuple) = pairs(x)
named_pairs(x::Vector{<:Pair{Symbol}}) = x

@inline Base.getproperty(x::TreeVectorAxes, leaf::Symbol) = _get_leaf(x, Val(leaf))
