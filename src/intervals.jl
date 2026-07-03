# intervals.jl
#
# A generic, immutable, statically-built *centered interval tree* for
# point-stabbing queries ("return every interval that contains a point").
#
# The tree is built once from a fixed set of intervals and then queried many
# times, which is exactly the liftover access pattern. It is fully parametric
# over the coordinate type `T` and the payload type `V`, so it carries no
# knowledge of genomics — that separation is the whole point of putting it in
# its own file.
#
# Design notes
#   * Immutable nodes with `Union{Nothing, IntervalNode{T,V}}` children. Julia
#     union-splits these small unions, so traversal stays type stable and
#     branch-predictable without any runtime `isa` checks.
#   * A node stores only the intervals that straddle its `center`, sorted by
#     `start`. Everything strictly left/right of `center` recurses into the
#     corresponding child. Correctness is independent of the exact `center`
#     value, so we pick an overflow-safe midpoint.
#   * Stabbing uses half-open interval semantics `start <= pos < stop`, which
#     is what genome coordinates need (a region's `stop` is one past its last
#     base).

"""
    Interval{T,V}(start, stop, value)

A half-open interval `[start, stop)` over coordinate type `T` carrying a
payload `value::V`.
"""
struct Interval{T,V}
    start::T
    stop::T
    value::V
end

"""
    IntervalNode{T,V}

Internal node of an [`IntervalTree`](@ref). `intervals` holds the intervals
straddling `center` (or, for a leaf, the node's whole bucket), sorted ascending
by `start`. `left`/`right` are the subtrees, or `nothing`.
"""
struct IntervalNode{T,V}
    center::T
    intervals::Vector{Interval{T,V}}
    left::Union{Nothing,IntervalNode{T,V}}
    right::Union{Nothing,IntervalNode{T,V}}
end

"""
    IntervalTree{T,V}

Immutable centered interval tree. Construct from a vector of intervals with
`IntervalTree{T,V}(intervals)`, or build an empty tree with `IntervalTree{T,V}()`.
"""
struct IntervalTree{T,V}
    root::Union{Nothing,IntervalNode{T,V}}
end

# Recursion guards. A bucket at or below `_MINBUCKET` intervals, or a branch
# that reaches `_MAXDEPTH`, becomes a leaf that is scanned linearly.
const _MINBUCKET = 32
const _MAXDEPTH = 32

IntervalTree{T,V}() where {T,V} = IntervalTree{T,V}(nothing)

function IntervalTree{T,V}(intervals::Vector{Interval{T,V}}) where {T,V}
    return IntervalTree{T,V}(_build(intervals, _MAXDEPTH))
end

function _build(intervals::Vector{Interval{T,V}}, depth::Int) where {T,V}
    isempty(intervals) && return nothing

    if depth == 0 || length(intervals) <= _MINBUCKET
        sort!(intervals; by = iv -> iv.start)
        return IntervalNode{T,V}(zero(T), intervals, nothing, nothing)
    end

    min_start = minimum(iv -> iv.start, intervals)
    max_stop = maximum(iv -> iv.stop, intervals)
    center = min_start + (max_stop - min_start) ÷ 2  # overflow-safe midpoint

    lefts = Interval{T,V}[]
    rights = Interval{T,V}[]
    here = Interval{T,V}[]
    for iv in intervals
        if iv.stop < center
            push!(lefts, iv)
        elseif iv.start > center
            push!(rights, iv)
        else
            push!(here, iv)
        end
    end
    sort!(here; by = iv -> iv.start)

    return IntervalNode{T,V}(center, here, _build(lefts, depth - 1), _build(rights, depth - 1))
end

"""
    stab(tree, pos) -> Vector{Interval{T,V}}

Return every interval containing `pos` under half-open semantics
`start <= pos < stop`. Allocates a fresh result vector per call.
"""
function stab(tree::IntervalTree{T,V}, pos::T) where {T,V}
    out = Interval{T,V}[]
    root = tree.root
    root === nothing && return out
    return _stab!(out, root, pos)
end

function _stab!(out::Vector{Interval{T,V}}, node::IntervalNode{T,V}, pos::T) where {T,V}
    ivals = node.intervals
    # `ivals` is sorted by start, so `first(ivals).start` is the minimum start;
    # if `pos` is below it, no interval in this bucket can contain `pos`.
    if !isempty(ivals) && !(pos < first(ivals).start)
        @inbounds for iv in ivals
            if iv.start <= pos < iv.stop
                push!(out, iv)
            end
        end
    end

    left = node.left
    if left !== nothing && pos <= node.center
        _stab!(out, left, pos)
    end

    right = node.right
    if right !== nothing && pos >= node.center
        _stab!(out, right, pos)
    end

    return out
end
