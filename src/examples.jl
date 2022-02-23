export examplegrv
import GenomicFeatures
import DataFrames

"""
    examplegrv(n::Int)

Mock up an example `GenomicRangeVector` with `n` entries.

```jldoctest
julia> using GenomicFeatureVectors

julia> x = examplegrv(10);

julia> length(x)
10

julia> seqinfo(x)[!,"seqname"]
3-element Vector{String}:
 "chrA"
 "chrB"
 "chrC"
```
"""
function examplegrv(n::Int)
    seqname = ["chrA", "chrB", "chrC"]
    len = [100, 200, 300]

    seqinfo = DataFrames.DataFrame(
        seqname = seqname,
        length = len,
        circular = [false, false, false],
        genome = ["mm10", "hg19", "rnor6"]
    )

    test = Vector{GenomicFeatures.Interval{Nothing}}(undef, n)
    for i in 1:n
        chosen = Integer(ceil(rand() * 3))
        l = len[chosen]
        start = Integer(ceil(rand() * l))
        last = min(start + Integer(ceil(rand() * 20)), l)
        test[i] = GenomicFeatures.Interval(seqname[chosen], start, last)
    end

    name = Vector{String}(undef, n)
    type = Vector{String}(undef, n)
    possible_types = ["protein_coding", "rRNA", "tRNA", "pseudogene", "lncRNA"]
    for i in 1:n
        name[i] = "Gene" * string(i)
        type[i] = possible_types[Integer(ceil(rand() * length(possible_types)))]
    end

    edata = DataFrames.DataFrame(name = name, type = type)
    meta = Dict("foo" => 2, "bar"=>"BAR");
    GenomicRangeVector(test, seqinfo, edata, meta)
end
