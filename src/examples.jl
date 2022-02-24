export examplegrv, exampleggrv
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

    edata = DataFrames.DataFrame("name" => name, "type" => type)
    meta = Dict("foo" => 2, "bar"=>"BAR");
    return GenomicRangeVector(test, seqinfo, edata, meta)
end

"""
    exampleggrv(ngroups::Int, nranges::Int)

Mock up an example `GroupedGenomicRangeVector` with `nranges` ranges split across `ngroups` groups.

```jldoctest
julia> using GenomicFeatureVectors

julia> x = exampleggrv(10, 20);

julia> length(x)
10

julia> seqinfo(x)[!,"seqname"]
3-element Vector{String}:
 "chrA"
 "chrB"
 "chrC"
```
"""
function exampleggrv(ngroups::Int, nranges::Int)
    rr = examplegrv(nranges)
    ed = elementdata(rr)
    ed = ed[!,[1]]

    # Make up some lengths.
    lengths = Vector{Int}(undef, ngroups)
    fill!(lengths, 0)
    for i in 1:nranges
        lengths[Integer(ceil(ngroups * rand()))] += 1
    end

    newnames = String[]
    for i in 1:ngroups
        for j in 1:lengths[i]
            push!(newnames, "Exon" * string(j))
        end
    end
    ed[!,"name"] = newnames
    setelementdata!(rr, ed)

    # Making up some metadata.
    name = Vector{String}(undef, ngroups)
    type = Vector{String}(undef, ngroups)
    possible_types = ["protein_coding", "rRNA", "tRNA", "pseudogene", "lncRNA"]
    for i in 1:ngroups
        name[i] = "Gene" * string(i)
        type[i] = possible_types[Integer(ceil(rand() * length(possible_types)))]
    end

    edata = DataFrames.DataFrame("name" => name, "type" => type)
    meta = Dict("foo" => 2, "bar" => "BAR");
    return GroupedGenomicRangeVector(rr, lengths, edata, meta; islength = true)
end
