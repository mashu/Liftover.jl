"""
    Liftover

Convert point coordinates between genome assemblies using UCSC chain files.

Motivated by [jeremymcrae/liftover](https://github.com/jeremymcrae/liftover)
(which traces back to pyliftover). A Julia rewrite — same liftover semantics
and a similar API, written from scratch with Julia-specific design. Chain files
are parsed into one centered interval tree per target contig for fast
strand-aware point queries.

```julia
using Liftover

converter = get_lifter("hg19", "hg38"; one_based=true)
converter["1"][103786442]                 # dictionary-style lift
convert_coordinate(converter, "1", 103786442)
query(converter, "1", 103786442)

# straight from a chain file
converter = ChainFile("/path/to/hg19ToHg38.over.chain.gz"; one_based=true)
converter["chr1"][103786442]
```
"""
module Liftover

using CodecZlib: GzipDecompressorStream
using Downloads: Downloads

export get_lifter, LiftOver, ChainFile, Target, Match, AbstractLifter
export query, convert_coordinate

include("intervals.jl")
include("headers.jl")
include("chain.jl")
include("target.jl")
include("chainfile.jl")
include("download.jl")
include("lifter.jl")

end # module
