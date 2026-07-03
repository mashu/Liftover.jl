# Getting started

## Installation

From the Julia REPL:

```julia
using Pkg
Pkg.add(url="https://github.com/mashu/Liftover.jl")
```

For local development, activate the package directory and instantiate:

```julia
using Pkg
Pkg.activate("path/to/Liftover.jl")
Pkg.instantiate()
```

## Creating a lifter

The simplest entry point is [`get_lifter`](@ref), which downloads a UCSC chain
file on first use and caches it under `~/.liftover`:

```julia
using Liftover

converter = get_lifter("hg19", "hg38")
```

[`LiftOver`](@ref) is an alias with the same behaviour (matching pyliftover):

```julia
lifter = LiftOver("hg19", "hg38")
```

To open a chain file you already have:

```julia
converter = ChainFile("/path/to/hg19ToHg38.over.chain.gz")
# or pass the path as the sole argument to get_lifter
converter = get_lifter("/path/to/hg19ToHg38.over.chain.gz")
```

## Querying coordinates

Three equivalent styles are supported:

```julia
converter["chr1"][pos]                       # dictionary-style
query(converter, "chr1", pos)                # explicit query
convert_coordinate(converter, "chr1", pos) # pyliftover-style
```

Each call returns a `Vector{Match}`. A [`Match`](@ref) holds the lifted contig,
position, and strand on the query genome:

```julia
matches = converter["chr1"][1000]
isempty(matches) || begin
    m = matches[1]
    println(m.contig, ' ', m.pos, ' ', m.strand)   # field access

    contig, pos, strand = m                          # tuple unpacking
    println(contig, ' ', pos, ' ', strand)
end
```

## Options

| Keyword          | Default                              | Description                                      |
| ---------------- | ------------------------------------ | ------------------------------------------------ |
| `one_based`      | `false`                              | Use 1-based coordinates for input and output   |
| `cache`          | `~/.liftover`                        | Directory for downloaded chain files             |
| `chain_server`   | `https://hgdownload.soe.ucsc.edu`    | UCSC mirror base URL (must preserve UCSC layout) |

```julia
converter = get_lifter("hg19", "hg38";
    one_based = true,
    cache = "/data/chain_cache",
    chain_server = "https://www.example.org",
)
```

The download URL follows UCSC's layout:

```
{chain_server}/goldenpath/{target}/liftOver/{target}To{Query}.over.chain.gz
```

## Multiple matches

Some positions map to more than one location (e.g. segmental duplications). The
result vector contains every match; an empty vector means the coordinate does
not lift over.

```julia
matches = converter["chr1"][pos]
length(matches)   # 0, 1, or more
```
