# chain.jl
#
# A single chain describes a run of aligned blocks between a target region and
# a query region. `ChainBuilder` consumes the header plus the alignment data
# lines that follow it, emitting one [`Interval`](@ref) per aligned block whose
# payload records where that block maps to in the query genome.
#
# Alignment data lines are either `size targetGap queryGap` (three fields) or a
# lone `size` on the final line of a chain.

"""
    Mapped

Payload stored on every target-space interval: where the block maps to in the
query genome (`start`/`stop`), the query contig id, the query strand
(`fwd_strand`), and the query contig size (needed to flip reverse-strand hits).
"""
struct Mapped
    start::Int64
    stop::Int64
    query_id::String
    fwd_strand::Bool
    size::Int64
end

"""
    parse_alignment(line) -> (size, target_gap, query_gap)

Parse a chain alignment data line. A single field means the final block, with
zero gaps. Three or more fields give `size`, `target_gap`, `query_gap` (extra
trailing fields are ignored, matching the reference implementation). Anything
else, or a non-integer field, throws `ArgumentError`.
"""
function parse_alignment(line::AbstractString)
    tokens = split(line)
    n = length(tokens)

    if n == 1
        size = tryparse(Int64, tokens[1])
        size === nothing && throw(ArgumentError("invalid alignment line: $line"))
        return size, Int64(0), Int64(0)
    elseif n >= 3
        size = tryparse(Int64, tokens[1])
        target_gap = tryparse(Int64, tokens[2])
        query_gap = tryparse(Int64, tokens[3])
        (size === nothing || target_gap === nothing || query_gap === nothing) &&
            throw(ArgumentError("invalid alignment line: $line"))
        return size, target_gap, query_gap
    else
        throw(ArgumentError("invalid alignment line: $line"))
    end
end

"""
    ChainBuilder

Mutable accumulator for one chain. Immutable header-derived fields are marked
`const`; only the running target/query cursors and the interval list mutate as
alignment lines are added.
"""
mutable struct ChainBuilder
    const target_id::String
    const query_id::String
    const fwd_strand::Bool
    const query_size::Int64
    const target_end::Int64
    const query_end::Int64
    target::Int64
    query::Int64
    const intervals::Vector{Interval{Int64,Mapped}}
end

function ChainBuilder(header::ChainHeader)
    return ChainBuilder(
        header.target_id,
        header.query_id,
        header.query_strand == '+',
        header.query_size,
        header.target_end,
        header.query_end,
        header.target_start,
        header.query_start,
        Interval{Int64,Mapped}[],
    )
end

ChainBuilder(header_line::AbstractString) = ChainBuilder(process_header(header_line))

"""
    add_line!(chain, line)

Parse one alignment data line and append the corresponding target-space
interval, then advance the target and query cursors past the block and its gaps.
"""
function add_line!(chain::ChainBuilder, line::AbstractString)
    size, target_gap, query_gap = parse_alignment(line)

    data = Mapped(chain.query, chain.query + size, chain.query_id, chain.fwd_strand, chain.query_size)
    push!(chain.intervals, Interval(chain.target, chain.target + size, data))

    chain.target += size + target_gap
    chain.query += size + query_gap
    return chain
end

"""
    validate(chain)

Check that the accumulated cursors landed exactly on the header's declared
target/query ends. Throws `ArgumentError` on mismatch.
"""
function validate(chain::ChainBuilder)
    chain.target == chain.target_end || throw(ArgumentError(
        "target end does not match expectations: $(chain.target) != $(chain.target_end)"))
    chain.query == chain.query_end || throw(ArgumentError(
        "query end does not match expectations: $(chain.query) != $(chain.query_end)"))
    return chain
end
