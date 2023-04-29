mutable struct TreeVectorAxes{NAMES, LEAVES}
    size::Int
    index::UnitRange{Int}
    leaves::LEAVES
end

function TreeVectorAxes(tree, idx = 1)
    names = Vector{Symbol}()
    size = 0
    leaves = Vector{Any}()
    firstidx = idx
    for (name, leaf) in named_pairs(tree)
        push!(names, name)
        if leaf isa Number
            throw(ArgumentError("Leaves of a tree vector should be vectors"))
        elseif leaf isa AbstractVector
            len = length(leaf)
            push!(leaves, TreeVectorAxes{(),Tuple{}}(len, idx:(idx+len-1), ()))
        else
            treeleaf = TreeVectorAxes(leaf, idx)
            len = length(treeleaf)
            push!(leaves, treeleaf)
        end
        size += len
        idx += len
    end
    _NAME = (names...,)
    _leaves = (leaves...,)
    _LEAF = typeof(_leaves)
    TreeVectorAxes{_NAME, _LEAF}(size, firstidx:(idx - 1), _leaves)
end

Base.isempty(::TreeVectorAxes) = false
Base.isempty(::TreeVectorAxes{0}) = true
Base.length(ax::TreeVectorAxes) = getfield(ax, :size)
Base.keys(::TreeVectorAxes{NAMES}) where NAMES = NAMES

named_pairs(ax::NamedTuple) = pairs(ax)
named_pairs(ax::Vector{<:Pair{Symbol}}) = ax
# TODO: named_pairs for TreeVectors

@inline Base.getproperty(ax::TreeVectorAxes, leaf::Symbol) = _get_leaf(ax, Val(leaf))

function _gen_get_leaf(::Type{TreeVectorAxes{NAME, LEAVES}}, ::Type{Val{LEAF}}) where {NAME, LEAVES, LEAF}
    idx = findfirst(==(LEAF), NAME)
    if idx === nothing
        return :(error("Unknown field: ", $LEAF))
    else
        return :(getindex(getfield(ax, :leaves), $idx))
    end
end

@generated function _get_leaf(ax::TreeVectorAxes, leaf)
    return _gen_get_leaf(ax, leaf)
end

function TreeVectorType(tree)
    T = Union{}
    for (name, leaf) in named_pairs(tree)
        if leaf isa AbstractVector
            T = promote_type(T, eltype(leaf))
        else
            T = promote_type(T, TreeVectorType(leaf))
        end
    end
    return T
end

struct TreeVector{T, V, A <: TreeVectorAxes} <: AbstractVector{T}
    data::V
    axes::A
end

function TreeVector(tree)
    axes = TreeVectorAxes(tree)
    T = TreeVectorType(tree)
    data = Vector{T}(undef, length(axes))
    tv = TreeVector{T, Vector{T}, typeof(axes)}(data, axes)
    copyto!(tv, tree)
end

function Base.copyto!(dest::TreeVector, tree)
    for (name, leaf) in named_pairs(tree)
        copyto!(getproperty(dest, name), leaf)
    end
    return dest
end

function Base.copyto!(dest::TreeVector, src::AbstractVector)
    copyto!(vec(dest), src)
    return dest
end

function Base.copyto!(::AbstractVector, ::Number)
    throw(ArgumentError("Cannot copy a number to a vector"))
end

function Base.copyto!(dest::TreeVector, bc::Base.Broadcast.Broadcasted)
    copyto!(vec(dest), bc)
end

function Base.getproperty(A::TreeVector, leaf::Symbol)
    ax = getproperty(getfield(A, :axes), leaf)
    return TreeVector{eltype(A), Vector{eltype(A)}, typeof(ax)}(getfield(A, :data), ax)
end

function Base.setproperty!(::TreeVector, ::Any, ::Symbol)
    throw(ArgumentError("Cannot change the structure of a TreeVector"))
end

Base.vec(A::TreeVector) = view(getfield(A, :data), getfield(getfield(A, :axes), :index))
Base.size(A::TreeVector) = (length(A),)
Base.length(A::TreeVector) = length(getfield(A, :axes))
Base.eltype(::TreeVector{T}) where {T} = T
Base.getindex(A::TreeVector, args...) = getindex(vec(A), args...)
Base.setindex!(A::TreeVector, v, args...) = setindex!(vec(A), v, args...)
Base.firstindex(A::TreeVector) = firstindex(vec(A))
Base.lastindex(A::TreeVector) = lastindex(vec(A))
Base.IndexStyle(::Type{<:TreeVector{T, V}}) where {T, V} = IndexStyle(V)
Base.iterate(A::TreeVector, state = firstindex(A)) = iterate(vec(A), state)

function Base.similar(A::TreeVector, ::Type{S}, dims::Dims)
    if dims == size(A)
        ax = getfield(A, :axes)
        return TreeVector{S, Vector{S}, typeof(ax)}(Vector{S}(undef, length(A)), ax)
    else
        throw(ArgumentError("Cannot change the structure of a TreeVector"))
    end
end

Base.keys(A::TreeVector) = keys(getfield(A, :axes))

# TODO: Implement Base.pairs for TreeVectors

TreeVectorAxes(A::TreeVector) = getfield(A, :axes)

# TODO: Implement arithmetic that preserves the tree structure
