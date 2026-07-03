using Test
using CodecZlib: GzipCompressorStream
using Liftover
using Liftover: process_header, parse_alignment, sanitize_prefix,
    Interval, IntervalTree, stab, Mapped

# Write `lines` to a temporary gzipped chain file and return its path.
function write_chain(lines::Vector{<:AbstractString})
    path = tempname() * ".chain.gz"
    open(GzipCompressorStream, path, "w") do io
        for line in lines
            write(io, line)
        end
    end
    return path
end

const MINIMAL = ["chain 0 chr1 10 + 0 10 chr1 10 + 10 30 2\n", "5 0 5\n", "5 0 5\n", "\n"]
# single reverse-strand block: target [1000,1100) -> query start 5000 on '-', qSize 8000
const REVERSE = ["chain 500 chr2 200000 + 1000 1100 chrX 8000 - 5000 5100 7\n", "100\n", "\n"]

@testset "Liftover.jl" begin

    @testset "header parsing" begin
        h = process_header("chain 0 chr1 10 + 0 10 chr2 20 - 3 13 5")
        @test h.target_id == "chr1"
        @test h.target_start == 0
        @test h.target_end == 10
        @test h.query_id == "chr2"
        @test h.query_size == 20
        @test h.query_strand == '-'
        @test h.query_start == 3
        @test h.query_end == 13

        @test_throws ArgumentError process_header("chain 0 chr1 10 + 0 10 chr1 10 + 10 30")      # 12 fields
        @test_throws ArgumentError process_header("notchain 0 chr1 10 + 0 10 chr1 10 + 10 30 2") # bad tag
        @test_throws ArgumentError process_header("chain 0 chr1 10 - 0 10 chr1 10 + 10 30 2")    # target strand
        @test_throws ArgumentError process_header("chain 0 chr1 10 + 0 10 chr1 10 ? 10 30 2")    # query strand
    end

    @testset "alignment parsing" begin
        @test parse_alignment("5000") == (5000, 0, 0)
        @test parse_alignment("5 0 5") == (5, 0, 5)
        @test parse_alignment("5 0 5 5") == (5, 0, 5)          # extra field ignored
        @test_throws ArgumentError parse_alignment("166661 50000")   # two fields
        @test_throws ArgumentError parse_alignment("s619 137 0")     # text in number
        @test_throws ArgumentError parse_alignment("invalid content")
    end

    @testset "interval tree stabbing" begin
        ivals = [Interval(0, 5, :a), Interval(5, 10, :b), Interval(2, 8, :c)]
        tree = IntervalTree{Int,Symbol}(copy(ivals))
        @test Set(iv.value for iv in stab(tree, 3)) == Set([:a, :c])   # 0<=3<5 and 2<=3<8
        @test Set(iv.value for iv in stab(tree, 5)) == Set([:b, :c])   # half-open: 5 excluded from :a
        @test isempty(stab(tree, 10))                                  # exclusive stop
        @test isempty(stab(IntervalTree{Int,Symbol}(), 3))             # empty tree
    end

    @testset "forward-strand lifting" begin
        cf = ChainFile(write_chain(MINIMAL))
        @test cf["chr1"][0] == [Match("chr1", 10, '+')]
        @test cf["chr1"][5] == [Match("chr1", 20, '+')]
        @test cf["chr1"][6] == [Match("chr1", 21, '+')]
        @test cf["chr1"][10] == Match[]                # boundary, exclusive stop
        @test cf["chr1"][6][1][2] == 21                # tuple-style indexing
    end

    @testset "reverse-strand lifting" begin
        cf = ChainFile(write_chain(REVERSE))
        @test cf["chr2"][1000] == [Match("chrX", 2999, '-')]
        @test cf["chr2"][1050] == [Match("chrX", 2949, '-')]
        @test cf["chr2"][1099] == [Match("chrX", 2900, '-')]
        @test cf["chr2"][1100] == Match[]              # boundary
    end

    @testset "one-based coordinates" begin
        fwd = ChainFile(write_chain(MINIMAL); one_based = true)
        @test fwd["chr1"][7] == [Match("chr1", 22, '+')]
        rev = ChainFile(write_chain(REVERSE); one_based = true)
        @test rev["chr2"][1051] == [Match("chrX", 2950, '-')]
    end

    @testset "prefix handling and missing contigs" begin
        cf = ChainFile(write_chain(REVERSE))          # file uses 'chr' prefix
        @test cf["2"][1050] == [Match("chrX", 2949, '-')]   # queried without prefix
        @test cf["chr2"].contig == cf["2"].contig == "chr2" # both resolve to chr2
        @test cf["chrZZ"][1050] == Match[]                  # missing contig -> empty
        @test "chr2" in keys(cf)
        @test !haskey(cf, "chrZZ")
    end

    @testset "Match tuple behaviour" begin
        m = Match("chr1", 21, '+')
        c, p, s = m
        @test (c, p, s) == ("chr1", 21, '+')
        @test m[1] == "chr1" && m[2] == 21 && m[3] == '+'
        @test length(m) == 3
        @test m.contig == "chr1" && m.pos == 21 && m.strand == '+'
    end

    @testset "convert_coordinate and query APIs" begin
        cf = ChainFile(write_chain(MINIMAL))
        @test query(cf, "chr1", 6) == [Match("chr1", 21, '+')]
        @test convert_coordinate(cf, "chr1", 6) == [Match("chr1", 21, '+')]
    end

    @testset "invalid and edge-case files" begin
        # target end does not match the header
        bad_end = ["chain 21270171362 chr1 249250621 + 10000 249233096 chr1 247249719 + 0 247199719 2\n",
                   "619 137 0\n", "166661 50000 50000\n", "\n"]
        @test_throws ArgumentError ChainFile(write_chain(bad_end))

        # truncated alignment line
        short = ["chain 21270171362 chr1 249250621 + 10000 249233096 chr1 247249719 + 0 247199719 2\n",
                 "619 137 0\n", "166661 50000\n", "\n"]
        @test_throws ArgumentError ChainFile(write_chain(short))

        # stray non-header content
        @test_throws ArgumentError ChainFile(write_chain(["invalid content"]))

        # missing final newline is fine
        no_newline = ["chain 0 chr1 10 + 0 10 chr1 10 + 10 30 2\n", "5 0 5\n", "5 0 5\n"]
        @test ChainFile(write_chain(no_newline))["chr1"][6] == [Match("chr1", 21, '+')]

        # comment lines are skipped
        comments = ["chain 0 chr1 10 + 0 10 chr1 10 + 10 30 2\n", "#5 0 5\n", "5 0 5\n", "5 0 5\n"]
        @test ChainFile(write_chain(comments))["chr1"][6] == [Match("chr1", 21, '+')]

        # very large coordinates near typemax(Int64)
        large = typemax(Int64) - 49
        big = ["chain 0 chr1 $large + 0 $large chr1 $large + 10 $(large + 20) 2\n",
               "5 0 5\n", "$(large - 5) 0 5\n", "\n"]
        @test ChainFile(write_chain(big))["chr1"][large - 50] == [Match("chr1", large - 40 + 5, '+')]
    end

    @testset "get_lifter from chain path" begin
        path = write_chain(MINIMAL)
        cf = get_lifter(path)              # no query build -> treat target as chain path
        @test cf["chr1"][6] == [Match("chr1", 21, '+')]
        @test_throws ArgumentError get_lifter("hg19")   # not a .chain.gz path
    end

end
