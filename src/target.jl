# target.jl
#
# A `Target` bundles every chain for one target contig into a single interval
# tree and answers coordinate queries against it. This is the object you reach
# after indexing a `ChainFile` by contig: `chainfile[contig][pos]`.

"""
    Match(contig, pos, strand)

A single lifted coordinate: the query `contig`, the lifted `pos`, and the query
`strand` (`'+'` or `'-'`).

For convenience it also behaves like the 3-tuple `(contig, pos, strand)` — it
is indexable and iterable — so `c, p, s = match` and `match[2]` both work, in
addition to the type-stable field access `match.contig`, `match.pos`,
`match.strand`.
"""
struct Match
    contig::String
    pos::Int64
    strand::Char
end

Base.length(::Match) = 3
Base.eltype(::Type{Match}) = Any

function Base.getindex(match::Match, i::Integer)
    i == 1 && return match.contig
    i == 2 && return match.pos
    i == 3 && return match.strand
    throw(BoundsError(match, i))
end

function Base.iterate(match::Match, state::Int = 1)
    state > 3 && return nothing
    return match[state], state + 1
end

Base.show(io::IO, match::Match) =
    print(io, '(', repr(match.contig), ", ", match.pos, ", ", repr(match.strand), ')')

"""
    Target(contig, tree, one_based)

Query object for a single target contig. Prefer building one from a vector of
chains with `Target(chains, one_based)`.
"""
struct Target
    contig::String
    tree::IntervalTree{Int64,Mapped}
    one_based::Bool
end

function Target(chains::Vector{ChainBuilder}, one_based::Bool)
    contig = chains[1].target_id

    total = 0
    for chain in chains
        total += length(chain.intervals)
    end

    intervals = Interval{Int64,Mapped}[]
    sizehint!(intervals, total)
    for chain in chains
        chain.target_id == contig || throw(ArgumentError("target ID mismatch"))
        append!(intervals, chain.intervals)
    end

    return Target(contig, IntervalTree{Int64,Mapped}(intervals), one_based)
end

"""
    query(target, pos) -> Vector{Match}

Lift coordinate `pos` through this contig's chains. Returns every match (there
can be zero, one, or several); the empty vector means the position does not lift.
"""
function query(target::Target, pos::Integer)
    shift = target.one_based ? Int64(1) : Int64(0)
    position = Int64(pos) - shift

    matches = Match[]
    root = target.tree.root
    root === nothing && return matches

    for interval in stab(target.tree, position)
        mapped = interval.value
        offset = position - interval.start
        remapped = mapped.start + offset
        mapped.fwd_strand || (remapped = mapped.size - remapped - 1)
        remapped += shift
        push!(matches, Match(mapped.query_id, remapped, mapped.fwd_strand ? '+' : '-'))
    end
    return matches
end

Base.getindex(target::Target, pos::Integer) = query(target, pos)
Base.show(io::IO, target::Target) = print(io, "Target(", repr(target.contig), ")")
