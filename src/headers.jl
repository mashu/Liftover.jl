# headers.jl
#
# Parsing and validation of a single chain header line. A UCSC chain header is
# 13 whitespace-separated fields:
#
#   chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
#
# We keep only the fields the lifter actually needs downstream.

"""
    ChainHeader

Parsed, validated fields of a chain header line that the lifter needs: the
target contig and its start/end, and the query contig, size, strand, and
start/end.
"""
struct ChainHeader
    target_id::String
    target_start::Int64
    target_end::Int64
    query_id::String
    query_size::Int64
    query_strand::Char
    query_start::Int64
    query_end::Int64
end

@inline function _header_int(field::AbstractString, line::AbstractString)
    value = tryparse(Int64, field)
    value === nothing && throw(ArgumentError("invalid header line: $line"))
    return value
end

"""
    process_header(line) -> ChainHeader

Split, validate, and parse a chain header line. Throws `ArgumentError` if the
line does not have 13 fields, does not start with `chain`, or has an invalid
strand.
"""
function process_header(line::AbstractString)
    fields = split(line)
    length(fields) == 13 || throw(ArgumentError("invalid header line: $line"))
    fields[1] == "chain" || throw(ArgumentError("header line does not start with 'chain': $line"))
    fields[5] == "+" || throw(ArgumentError("target strand is not '+': $line"))
    (fields[10] == "+" || fields[10] == "-") ||
        throw(ArgumentError("query strand is not '+' or '-': $line"))

    return ChainHeader(
        String(fields[3]),
        _header_int(fields[6], line),
        _header_int(fields[7], line),
        String(fields[8]),
        _header_int(fields[9], line),
        first(fields[10]),
        _header_int(fields[11], line),
        _header_int(fields[12], line),
    )
end
