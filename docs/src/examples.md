# Examples

These examples show typical workflows with genomic coordinates. They use
`one_based=true` because VCF and most tabular variant data are 1-based.

## Lift a single variant position (hg19 → hg38)

```julia
using Liftover

hg19_to_hg38 = get_lifter("hg19", "hg38"; one_based=true)

chrom = "1"
pos = 103786442   # 1-based position on hg19

matches = hg19_to_hg38[chrom][pos]

if isempty(matches)
    println("Position does not lift over")
else
    for m in matches
        println("hg38: ", m.contig, ":", m.pos, " (strand ", m.strand, ")")
    end
end
```

## Batch-lift a table of SNP positions

```julia
using Liftover

converter = get_lifter("hg19", "hg38"; one_based=true)

snps = [
    ("1", 103786442),
    ("2", 234567890),
    ("X", 154331033),
]

for (chrom, pos) in snps
    matches = converter[chrom][pos]
    if isempty(matches)
        println(chrom, ":", pos, " -> (no liftover)")
    else
        m = matches[1]   # take first hit when unambiguous
        println(chrom, ":", pos, " -> ", m.contig, ":", m.pos)
    end
end
```

When a position has multiple mappings, inspect the full vector or apply your
own disambiguation rule (e.g. keep the match on the same contig).

## Lift coordinates from a local chain file

Useful when working offline or with custom chain files:

```julia
using Liftover

chain_path = "/data/liftover/mm10ToMm39.over.chain.gz"
converter = ChainFile(chain_path; one_based=true)

# mouse chromosome 1, 1-based
matches = converter["chr1"][12345678]
```

## Round-trip check (hg38 → hg19 → hg38)

Verify that a position maps consistently through a round trip:

```julia
using Liftover

to_hg38 = get_lifter("hg19", "hg38"; one_based=true)
to_hg19 = get_lifter("hg38", "hg19"; one_based=true)

chrom, pos = "1", 103786442

fwd = to_hg38[chrom][pos]
@assert !isempty(fwd)
m = fwd[1]

back = to_hg19[m.contig][m.pos]
@assert !isempty(back)
b = back[1]

println("Original:  chr", chrom, ":", pos)
println("Via hg38:  ", m.contig, ":", m.pos)
println("Back:      ", b.contig, ":", b.pos)
```

## Use a UCSC mirror

If the default UCSC download host is slow or blocked, point at a mirror that
mirrors the same directory layout:

```julia
using Liftover

converter = get_lifter("hg19", "hg38";
    chain_server = "https://hgdownload.soe.ucsc.edu",
    cache = joinpath(homedir(), ".liftover"),
)
```

## 0-based BED-style coordinates

For BED files (0-based start, half-open end), use the default coordinate system:

```julia
using Liftover

converter = get_lifter("hg19", "hg38")   # one_based=false by default

# BED record: chr1  1000  2000  (covers bases at indices 1000..1999)
bed_start = 1000
bed_end_exclusive = 2000

start_lift = converter["chr1"][bed_start]
end_lift = converter["chr1"][bed_end_exclusive - 1]   # last included base

for (label, matches) in [("start", start_lift), ("end", end_lift)]
    if isempty(matches)
        println(label, ": no liftover")
    else
        m = matches[1]
        println(label, ": ", m.contig, ":", m.pos)
    end
end
```

Lift start and end independently; do not assume the lifted interval length is
preserved across assemblies.

## Inspect available contigs

```julia
using Liftover

converter = get_lifter("hg19", "hg38")

for contig in sort(collect(keys(converter)))
    println(contig)
end
```
