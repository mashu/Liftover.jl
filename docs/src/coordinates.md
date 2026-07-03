# Coordinates

Understanding coordinate conventions is essential when lifting genomic positions
between assemblies.

## 0-based vs 1-based

By default Liftover.jl uses **0-based, half-open** coordinates — the UCSC and BED
convention where the first base on a contig is position `0` and intervals are
`[start, stop)`.

Pass `one_based=true` when creating the lifter to work in **1-based** coordinates
(VCF, GFF, SAM POS field, and most human-facing reports):

```julia
# 0-based (default) — BED-style
converter = ChainFile("hg19ToHg38.over.chain.gz")
converter["chr1"][0]    # first base

# 1-based — VCF-style
converter = ChainFile("hg19ToHg38.over.chain.gz"; one_based=true)
converter["chr1"][1]    # first base
```

The `one_based` flag applies symmetrically: you query in the chosen system and
receive results in the same system.

## Contig naming (`chr` prefix)

UCSC chain files typically use `chr`-prefixed contig names (`chr1`, `chrX`).
Other pipelines omit the prefix (`1`, `X`). Liftover reconciles these
automatically:

```julia
converter = get_lifter("hg19", "hg38")

converter["chr1"][pos]   # works when the chain file uses chr1
converter["1"][pos]      # same result — prefix added or stripped as needed
```

If the contig is not present in the chain file, queries return an empty vector
rather than raising an error:

```julia
converter["chrZZ"][1000] == Match[]
```

Use `keys` and `haskey` on the lifter to inspect available contigs:

```julia
"chr1" in keys(converter)
haskey(converter, "chr1")
```

## Strand

Each [`Match`](@ref) reports the strand on the **query** (destination) genome:

- `'+'` — forward strand alignment
- `'-'` — reverse strand alignment (position is remapped accordingly)

```julia
m = converter["chr2"][1050][1]
m.strand   # '+' or '-'
```

## Interval boundaries

Coordinates use half-open intervals: a chain block covering target positions
`[1000, 1100)` lifts positions `1000` through `1099`. Position `1100` is
outside the block and returns no match unless covered by another chain.

This matches UCSC liftOver behaviour and is important when lifting BED intervals
— lift the start and end positions separately, or use a dedicated interval
liftover tool for full interval semantics.
