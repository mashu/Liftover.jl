# Liftover.jl

Convert point coordinates between genome assemblies using UCSC chain files.

This is a Julia port of [jeremymcrae/liftover](https://github.com/jeremymcrae/liftover)
(a fast Cython/C++ implementation inspired by pyliftover). Chain files are parsed
into one centered interval tree per target contig, giving fast, strand-aware
point queries with a dictionary-style interface.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/<you>/Liftover.jl")
```

## Usage

```julia
using Liftover

converter = get_lifter("hg19", "hg38"; one_based=true)
chrom = "1"
pos = 103786442

converter[chrom][pos]                       # dictionary-style
convert_coordinate(converter, chrom, pos)   # pyliftover-style
query(converter, chrom, pos)                # synonym

# straight from a chain file
converter = ChainFile("/path/to/hg19ToHg38.over.chain.gz"; one_based=true)
converter["chr1"][pos]

# use a UCSC mirror (must preserve the UCSC URL layout)
converter = get_lifter("hg19", "hg38"; chain_server="https://www.example.org")

# pyliftover-style entry-point name
lifter = LiftOver("hg19", "hg38")           # or LiftOver(path_to_chain)
```

Each query returns a `Vector{Match}`. A `Match` exposes type-stable fields and
also unpacks like the tuple `(contig, pos, strand)`:

```julia
matches = converter["1"][103786442]
for m in matches
    println(m.contig, ' ', m.pos, ' ', m.strand)
end

contig, pos, strand = matches[1]            # tuple-style destructuring
```

An empty vector means the position does not lift over. When the queried contig
uses a different `chr` prefix convention than the chain file, the prefix is
reconciled automatically (`"1"` finds `"chr1"` and vice versa).

## Coordinates

By default coordinates are 0-based (the UCSC/BED convention). Pass
`one_based=true` to query with, and receive back, 1-based coordinates.

## Design

The package follows one responsibility per file:

| File            | Responsibility                                                  |
| --------------- | --------------------------------------------------------------- |
| `intervals.jl`  | Generic, parametric centered interval tree and point stabbing.  |
| `headers.jl`    | Chain header parsing and validation.                            |
| `chain.jl`      | Alignment-line parsing and per-chain interval accumulation.     |
| `target.jl`     | Per-contig query object and strand-aware coordinate remapping.  |
| `chainfile.jl`  | Chain-file reading, prefix handling, and the `ChainFile` lifter.|
| `download.jl`   | Atomic download of chain files.                                 |
| `lifter.jl`     | `get_lifter` / `LiftOver` convenience constructors.             |

The interval tree is parametric over coordinate and payload types (`IntervalTree{T,V}`),
built once and queried many times. Traversal relies on multiple dispatch and
small-`Union` splitting (`Union{Nothing, IntervalNode{T,V}}`) to stay type
stable without runtime type checks. Genomics-specific types (`Mapped`, `Match`)
are kept out of the tree so the tree is reusable on its own.

## Tests

```julia
using Pkg
Pkg.test("Liftover")
```

## License

MIT, matching the original project.
