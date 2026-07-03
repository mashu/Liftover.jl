# lifter.jl
#
# Convenience constructor mirroring the Python `get_lifter` / `LiftOver` entry
# points: name two genome builds and get back a ready `ChainFile`, downloading
# and caching the UCSC chain file on first use, or point it straight at a
# `.chain.gz` path.

"""
    get_lifter(target, query=nothing; cache=nothing, one_based=false,
               chain_server="https://hgdownload.soe.ucsc.edu") -> ChainFile

Create a converter mapping coordinates from the `target` genome build to the
`query` build (e.g. `get_lifter("hg19", "hg38")`). The chain file is downloaded
into `cache` (default `~/.liftover`) if not already present.

If `query` is `nothing`, `target` is treated as a path to a `.chain.gz` file and
opened directly.

`chain_server` lets you point at a UCSC mirror; it must preserve the UCSC URL
layout, i.e. `{chain_server}/goldenpath/{target}/liftOver/{target}To{Query}.over.chain.gz`.
"""
function get_lifter(target::AbstractString, query = nothing;
                    cache = nothing,
                    one_based::Bool = false,
                    chain_server::AbstractString = "https://hgdownload.soe.ucsc.edu")
    cache_dir = cache === nothing ? joinpath(homedir(), ".liftover") : String(cache)
    mkpath(cache_dir)

    if query === nothing
        endswith(target, ".chain.gz") ||
            throw(ArgumentError("target must be a chain file if no query is provided"))
        chain_path = String(target)
    else
        target_build = lowercasefirst(String(target))
        query_build = uppercasefirst(String(query))
        basename = string(target_build, "To", query_build, ".over.chain.gz")
        chain_path = joinpath(cache_dir, basename)
        if !isfile(chain_path)
            url = string(chain_server, "/goldenpath/", target_build, "/liftOver/", basename)
            download_chain(url, chain_path)
        end
    end

    return ChainFile(chain_path; one_based = one_based)
end

"""
    LiftOver(target, query=nothing; kwargs...)

Alias for [`get_lifter`](@ref), matching the pyliftover entry-point name.
"""
const LiftOver = get_lifter
