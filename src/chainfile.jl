# chainfile.jl
#
# Parses a (optionally gzipped) UCSC chain file into one interval tree per
# target contig, and exposes the pyliftover-style access surface:
#
#   chainfile[contig][pos]              # dictionary-style
#   query(chainfile, contig, pos)
#   convert_coordinate(chainfile, contig, pos)
#
# `AbstractLifter` exists so the generic `convert_coordinate` fallback is shared
# by dispatch rather than duplicated; `ChainFile` is the one concrete lifter here.

"""
    AbstractLifter

Supertype for coordinate lifters. A concrete lifter must implement
`query(lifter, contig, pos)`; `convert_coordinate` is then provided for free.
"""
abstract type AbstractLifter end

"""
    convert_coordinate(lifter, contig, pos) -> Vector{Match}

Lift `pos` on `contig` from the target genome to the query genome. Alias for
[`query`](@ref); provided for compatibility with pyliftover.
"""
convert_coordinate(lifter::AbstractLifter, contig::AbstractString, pos::Integer) =
    query(lifter, contig, pos)

# Open `path`, transparently decompressing if it carries the gzip magic bytes,
# so both `.chain.gz` and plain `.chain` files just work.
function _open_chain_stream(path::AbstractString)
    io = open(path)
    magic = read(io, 2)
    seekstart(io)
    if length(magic) == 2 && magic[1] == 0x1f && magic[2] == 0x8b
        return GzipDecompressorStream(io)
    end
    return io
end

"""
    open_chainfile(path, one_based) -> Dict{String,Target}

Read a chain file and build a `Target` per target contig. Chains are separated
by blank lines; `#` lines are comments; every other non-header line is alignment
data for the current chain. Each completed chain is validated before use.
"""
function open_chainfile(path::AbstractString, one_based::Bool)
    chains = Dict{String,Vector{ChainBuilder}}()
    current = nothing

    stream = _open_chain_stream(path)
    try
        for raw in eachline(stream)
            line = rstrip(raw, '\r')

            if isempty(line)
                if current !== nothing
                    validate(current)
                    push!(get!(() -> ChainBuilder[], chains, current.target_id), current)
                    current = nothing
                end
            elseif line[1] == '#'
                continue
            elseif startswith(line, "chain")
                current = ChainBuilder(line)
            else
                if current === nothing
                    # Reproduce the reference behaviour: a stray data line before
                    # any header is parsed as alignment data, which raises
                    # "invalid alignment line" for non-numeric content.
                    parse_alignment(line)
                    throw(ArgumentError("alignment line before any chain header: $line"))
                end
                add_line!(current, line)
            end
        end
    finally
        close(stream)
    end

    if current !== nothing
        validate(current)
        push!(get!(() -> ChainBuilder[], chains, current.target_id), current)
    end

    targets = Dict{String,Target}()
    for (contig, builders) in chains
        targets[contig] = Target(builders, one_based)
    end
    return targets
end

"""
    sanitize_prefix(contig, target_prefixed) -> String

Reconcile a `chr` prefix between a queried `contig` and the file's convention,
so `"1"` finds `"chr1"` (and vice versa) when only one side uses the prefix.
"""
function sanitize_prefix(contig::AbstractString, target_prefixed::Bool)
    contig_prefixed = startswith(contig, "chr")
    contig_prefixed == target_prefixed && return String(contig)
    return contig_prefixed ? String(contig[4:end]) : string("chr", contig)
end

"""
    ChainFile(path; one_based=false) -> ChainFile

Open a chain file for lifting coordinates from the target genome to the query
genome. Set `one_based=true` if you query with (and want back) 1-based
coordinates.
"""
struct ChainFile <: AbstractLifter
    path::String
    targets::Dict{String,Target}
    target_prefixed::Bool
    one_based::Bool
    empty_target::Target
end

function ChainFile(path::AbstractString; one_based::Bool = false)
    targets = open_chainfile(path, one_based)
    target_prefixed = any(startswith("chr"), keys(targets))
    empty_target = Target("", IntervalTree{Int64,Mapped}(), one_based)
    return ChainFile(String(path), targets, target_prefixed, one_based, empty_target)
end

function Base.getindex(chainfile::ChainFile, contig::AbstractString)
    haskey(chainfile.targets, contig) && return chainfile.targets[contig]
    key = sanitize_prefix(contig, chainfile.target_prefixed)
    return get(chainfile.targets, key, chainfile.empty_target)
end

"""
    query(chainfile, contig, pos) -> Vector{Match}

Lift `pos` on `contig` from the target genome to the query genome.
"""
query(chainfile::ChainFile, contig::AbstractString, pos::Integer) = chainfile[contig][pos]

Base.keys(chainfile::ChainFile) = keys(chainfile.targets)
Base.haskey(chainfile::ChainFile, contig::AbstractString) = haskey(chainfile.targets, contig)
Base.show(io::IO, chainfile::ChainFile) = print(io, "ChainFile(", repr(chainfile.path), ")")
