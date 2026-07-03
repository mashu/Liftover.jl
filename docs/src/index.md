# Liftover.jl

Liftover.jl converts point coordinates between genome assemblies using
[UCSC chain files](https://genome.ucsc.edu/goldenPath/help/chain.html).
Motivated by [jeremymcrae/liftover](https://github.com/jeremymcrae/liftover)
(which traces back to [pyliftover](https://github.com/jeremymcrae/pyliftover)):
a Julia rewrite of the same idea, written from scratch rather than translated
from the original code.

Chain files are parsed into one centered interval tree per target contig, giving
fast, strand-aware point queries with a dictionary-style interface.

## Quick example

```julia
using Liftover

converter = get_lifter("hg19", "hg38"; one_based=true)
matches = converter["1"][103786442]

for m in matches
    println(m.contig, ' ', m.pos, ' ', m.strand)
end
```

An empty result means the position does not lift over. See
[Getting started](@ref) and [Examples](@ref) for more detail.

## Features

- **Fast queries** — interval trees built once, queried many times.
- **Strand-aware** — reverse-strand alignments are remapped correctly.
- **Flexible input** — UCSC build names, local `.chain.gz` files, or custom mirrors.
- **Coordinate systems** — 0-based (UCSC/BED) or 1-based (VCF/GFF) via `one_based`.
- **Prefix handling** — `"1"` and `"chr1"` are reconciled automatically.

```@index
Pages = ["getting_started.md", "coordinates.md", "examples.md", "api.md"]
```
