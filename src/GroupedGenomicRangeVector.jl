export GroupedGenomicRangeVector
export intervals, setintervals!
export seqinfo, setseqinfo!
export ranges, setranges!
export lengths

############# Class definition ################

function create_cumulative(groups::AbstractVector{<:Integer}, islength::Bool)
    cumulative = Vector{Int}(undef, length(groups))

    if length(groups) > 0
        if islength
            cumulative[1] = groups[1]
            for i in 2:length(groups)
                if groups[i] < 0
                    throw(ErrorException("lengths in 'groups' cannot be negative"))
                end
                cumulative[i] = cumulative[i-1] + groups[i]
            end
        else
            for i in 1:length(groups)
                if i > 1 && groups[i] < groups[i-1]
                    throw(ErrorException("endpoints in 'groups' must be non-decreasing"))
                end
                cumulative[i] = groups[i]
            end
        end
    end

    return cumulative
end

mutable struct GroupedGenomicRangeVector{T} <: GenomicFeatureVector{GenomicRangeVector{T}}
    ranges::GenomicRangeVector{T}
    cumulative::Vector{Int}
    elementdata::DataFrames.DataFrame
    metadata::Dict{String,Any}

    function GroupedGenomicRangeVector(T)
        mockelementdata = mock_elementdata(0)
        new{T}(GenomicRangeVector(T), Int[], mockelementdata, Dict{String,Any})
    end

    function GroupedGenomicRangeVector(ranges::GenomicRangeVector{T}, groups::Vector{Int}; islength = true) where {T}
        cumulative = create_cumulative(groups, islength)
        mockelementdata = mock_elementdata(length(groups))
        new{T}(ranges, cumulative, mockelementdata, Dict{String,Any})
    end

    function GroupedGenomicRangeVector(
            ranges::GenomicRangeVector{T}, 
            groups::Vector{Int}, 
            elementdata::DataFrames.DataFrame, 
            metadata::Dict{String,Any}; 
            islength = true
        ) where {T}

        cumulative = create_cumulative(groups, islength)
        new{T}(ranges, cumulative, elementdata, metadata)
    end
end

############# Custom methods ################
#
"""
    intervals(x::GroupedGenomicRangeVector{T}) 

Return the underlying vector of `GenomicFeatures.Intervals`.
This is equivalent to the concatenation of the result of calling `intervals` on each group.

# Examples
```jldoctest
julia> using GenomicFeatureVectors

julia> x = exampleggrv(10, 20);

julia> length(x)
10

julia> iv = intervals(x);

julia> length(iv)
20
```
"""
function intervals(x::GroupedGenomicRangeVector{T}) where {T}
    return intervals(x.ranges)
end

"""
    setintervals!(x::GroupedGenomicRangeVector{T}, value::Vector{GenomicFeatures.Interval{T}}) 

Set the underlying vector of `GenomicFeatures.Intervals` to `value`.
`value` should contain the concatenation of intervals across all groups in `x`,
and should have the same length as `intervals(x)`.

# Examples
```jldoctest
julia> using GenomicFeatureVectors, GenomicFeatures

julia> x = exampleggrv(10, 20);

julia> vec = copy(intervals(x));

julia> for (i, v) in enumerate(vec)
           vec[i] = GenomicFeatures.Interval(
               replace(v.seqname, "chr" => "whee"), 
               v.first, 
               v.last, 
               v.strand, 
               v.metadata
           )
       end

julia> setintervals!(x, vec);
"""
function setintervals!(x::GroupedGenomicRangeVector{T}, intervals::GenomicRangeVector{T}) where {T}
    setintervals!(x, intervals)
    return 
end

"""
    seqinfo(x::GroupedGenomicRangeVector{T}) 

Return a `DataFrame` containing the sequence information for the ranges in `x`.
This is the same across all groups in `x`.

# Examples
```jldoctest
julia> using GenomicFeatureVectors

julia> x = exampleggrv(10, 20);

julia> names(seqinfo(x))
4-element Vector{String}:
 "seqname"
 "length"
 "circular"
 "genome"
```
"""    
function seqinfo(x::GroupedGenomicRangeVector{T}; check = true) where {T}
    return seqinfo(x.ranges)
end

"""
    setseqinfo!(x::GroupedGenomicRangeVector{T}, value::DataFrames.DataFrame) 

Set the sequence information for all groups in `x` to the `DataFrame` `value`.
This should contain the same columns as described in the `seqinfo` method for `GenomicRangeVector` objects.

# Examples
```jldoctest
julia> using GenomicFeatureVectors

julia> x = exampleggrv(10, 20);

julia> si = copy(seqinfo(x));

julia> si[!,"length"] .*= 10; # modifying the sequence lengths

julia> setseqinfo!(x, si);
```
"""
function setseqinfo!(x::GroupedGenomicRangeVector{T}, value::DataFrames.DataFrame) where {T}
    setseqinfo!(x.ranges, value)
    return x
end

"""
    ranges(x::GroupedGenomicRangeVector{T})

Return the underlying `GenomicRangeVector`,
containing the concatenation of ranges from all groups in `x`.

# Examples
```jldoctest
julia> using GenomicFeatureVectors

julia> x = exampleggrv(10, 20);

julia> y = ranges(x);

julia> length(y)
20
```
"""
function ranges(x::GroupedGenomicRangeVector{T}) where {T}
    return x.ranges
end

"""
    setranges!(x::GroupedGenomicRangeVector{T}, value::GenomicRangeVector{T})

Set the underlying `GenomicRangeVector` in `x` to `value`.
`value` should contain the concatenation of ranges across all groups in `x`,
and should have the same length as `ranges(x)`.

# Examples
```jldoctest
julia> using GenomicFeatureVectors

julia> x = exampleggrv(10, 20);

julia> y = copy(ranges(x));

julia> df = elementdata(y);

julia> df[!,"name"] = [replace(n, "Exon"=>"FOO") for n in df[!,"name"]];

julia> setelementdata!(y, df);

julia> setranges!(x, y);
```
"""
function setranges!(x::GroupedGenomicRangeVector{T}, value::GenomicRangeVector{T}) where {T}
    if length(value) != length(x.ranges)
        throw(ErrorException("lengths of 'value' and 'ranges(x)' must be the same"))
    end
    x.ranges = value
    return x
end

"""
    lengths(x::GroupedGenomicRangeVector{T})

Return a vector of length equal to `x`, containing the length of each group. 

# Examples
```jldoctest
julia> using GenomicFeatureVectors

julia> x = exampleggrv(10, 20);

julia> l = lengths(x);

julia> length(l)
10
```
"""
function lengths(x::GroupedGenomicRangeVector{T}) where {T}
    output = Vector{Int}(undef, length(x))
    if length(output) > 0
        output[1] = x.cumulative[1]
        for i in 2:length(output)
            output[i] = x.cumulative[i] - x.cumulative[i-1]
        end
    end
    return output
end

############# Miscellaneous methods ################

"""
    length(x::GroupedGenomicRangeVector{T})

Return the number of groups in `x`.
To get the size of each group, use [`lengths`](@ref) instead.
"""
function Base.length(x::GroupedGenomicRangeVector{T}) where {T}
    return length(x.cumulative)
end

function single_element_range(x::GroupedGenomicRangeVector{T}, i::Int) where {T}
    s = (i == 1 ? 0 : x.cumulative[1]) + 1
    l = x.cumulative[i]
    return s:l
end

"""
    getindex(x::GroupedGenomicRangeVector{T}, i::Int)

Return a `GenomicRangeVector{T}` for the group at index `i` in `x`.

# Examples
```jldoctest
julia> using GenomicFeatureVectors;

julia> x = exampleggrv(10, 20);

julia> y = x[1];

julia> y = x[2];

julia> typeof(y)
GenomicRangeVector{Nothing}
```
"""
function Base.getindex(x::GroupedGenomicRangeVector{T}, i::Int) where {T}
    rr = single_element_range(x, i)
    return x.ranges[rr]
end

"""
    setindex!(x::GroupedGenomicRangeVector{T}, value::GenomicRangeVector{T}, i::Int)

Set the group at index `i` to `value`, returning a reference to the modified object.
This is done efficiently if `value` and `x[i]` are of the same length.

# Examples
```jldoctest
julia> using GenomicFeatureVectors;

julia> x = exampleggrv(10, 20);

julia> x[1] = x[1]; # same length

julia> x[2] = x[1]; # possibly different length

julia> length(x[2]) == length(x[1])
true
```
"""
function Base.setindex!(x::GroupedGenomicRangeVector{T}, value::GenomicRangeVector{T}, i::Int) where {T}
    rr = single_element_range(x, i)
    if length(value) == length(rr)
        x.ranges[rr] = value
    else
        # Reassembly required if lengths don't match up.
        before = 1:(first(rr) - 1)
        after = (last(rr) + 1):length(x.ranges)
        x.ranges = vcat(x.ranges[before], value, x.ranges[after])
        delta = length(value) - length(rr)
        for j in i:length(x.cumulative)
            x.cumulative[j] += delta
        end
    end
    return x
end

# Defining a more efficient method when the requested indices are contiguous.
function getindex_range(x::GroupedGenomicRangeVector, I::AbstractRange{<:Integer})
    if length(I) > 0
        rr1 = first(single_element_range(x, first(I)))
        rr2 = last(single_element_range(x, last(I)))
        sub = x.cumulative[I] .- (rr1 - 1)
        return (rr1:rr2, sub)
    else
        return (1:0, Int[])
    end
end

function getindex_range(x::GroupedGenomicRangeVector, I::AbstractVector{<:Integer})
    new_range_index = Int[]
    latest = 0
    new_cum = Int[]
    for i in I
        rr = single_element_range(x, i)
        append!(new_range_index, rr)
        latest += length(rr)
        push!(new_cum, latest)
    end
    return (new_range_index, new_cum)
end

"""
    getindex(x::GroupedGenomicRangeVector{T}, I)

Take a slice of `x` at the specified indices `I`, returning a new `GroupedGenomicRangeVector`.
This can be done most efficiently if `I` is a range.

# Examples
```jldoctest
julia> using GenomicFeatureVectors;

julia> x = exampleggrv(10, 20);

julia> y = x[2:5];

julia> length(y)
4

julia> typeof(y)
GroupedGenomicRangeVector{Nothing}

julia> y = x[[3,4,5,1,2,3]];

julia> length(y)
6
```
"""
function Base.getindex(x::GroupedGenomicRangeVector{T}, I) where {T}
    idx, = to_indices(x, (I,))
    new_range_index, new_cum = getindex_range(x, idx)
    return GroupedGenomicRangeVector(x.ranges[new_range_index], new_cum, x.elementdata[idx,:], x.metadata; islength = false)
end

function setindex_internal(x::GroupedGenomicRangeVector{T}, value::GroupedGenomicRangeVector{T}, I::AbstractRange{<:Integer}) where {T}
    rr1 = first(single_element_range(x, first(I)))
    rr2 = last(single_element_range(x, last(I)))
    rr = rr1:rr2

    if length(rr) == length(value.ranges)
        # Making copies to avoid leaving the object in a broken
        # state if there are failures down the line.
        copyranges = copy(x.ranges)
        copyranges[rr] = value.ranges

        copycumulative = copy(x.cumulative)
        delta = first(rr) - 1
        for (i, j) in enumerate(I)
            copycumulative[j] = value.cumulative[i] + delta
        end

        return (copyranges, copycumulative)
    else
        # Some reassembly required.
        before = 1:(first(rr) - 1)
        after = (last(rr) + 1):length(x.ranges)
        copyranges = vcat(x.ranges[before], value.ranges, x.ranges[after])

        copycumulative = copy(x.cumulative)

        delta1 = first(rr) - 1
        delta2 = length(value.ranges) - length(rr)

        for (i, j) in enumerate(I)
            copycumulative[j] = value.cumulative[i] + delta1
        end

        for j in (last(I)+1):length(x.cumulative)
            copycumulative[j] += delta2
        end

        return (copyranges, copycumulative)
    end
end

function setindex_internal(x::GroupedGenomicRangeVector{T}, value::GroupedGenomicRangeVector{T}, I::AbstractVector{<:Integer}) where {T}
    # Building a back-mapping index.
    mapping = Dict{Int,Int}()
    for (i, j) in enumerate(I)
        mapping[j] = i
    end

    new_range = Vector{GenomicRangeVector{T}}(undef, length(x))
    new_cum = Vector{Int}(undef, length(x))
    sofar = 0

    for i in 1:length(x)
        local curset::GenomicRangeVector{T}
        if haskey(mapping, i)
            curset = value[mapping[i]]
        else
            curset = x[i]
        end

        new_range[i] = curset
        sofar += length(curset)
        new_cum[i] = sofar
    end

    return (vcat(new_range...), new_cum)
end

"""
    setindex!(x::GroupedGenomicRangeVector{T}, value::GroupedGenomicRangeVector{T}, I)

Set `x` to `value` at the specified indices `I`, returning a reference to the modified object.
This can be done most efficiently if `I` is a range.

# Examples
```jldoctest
julia> using GenomicFeatureVectors;

julia> x = exampleggrv(10, 20);

julia> x[2:5] = x[2:5];

julia> x[2:5] = x[4:7];

julia> x[[3,4,5,1,2,3]] = x[1:6];
```
"""
function Base.setindex!(x::GroupedGenomicRangeVector{T}, value::GroupedGenomicRangeVector{T}, I) where {T}
    idx, = to_indices(x, (I,))
    if length(idx) == 0
        return # no-op
    end

    new_range, new_cum = setindex_internal(x, value, idx)

    # If metadata replacement succeeds, we can replace the other fields
    # under the assumption that they are now no-throws.
    x.elementdata[idx,:] = value.elementdata

    x.ranges = new_range
    x.cumulative = new_cum
    return x
end

"""
    vcat(A::Vararg{GroupedGenomicRangeVector{T}})

Returns a `GroupedGenomicRangeVector{T}` containing the concatenated groupings from all `A`.
The length of the output is equal to the sum of the lengths of all input objects,
and the [`lengths`](@ref) of the output is equal to the concatenation of the [`lengths](@ref) of the inputs.

The underlying [`ranges`](@ref) are concatenated as described in the `vcat` method for `GenomicRangeVector` objects.

Element data is concatenated by row, so each entry of `A` should have the same type and name of columns in their [`elementdata`](@ref).

Metadata in [`metadata`](@ref) is combined across `A`.
Earlier occurrences of each key take precedence, and later occurrences are dropped silently.

# Examples
```jldoctest
julia> using GenomicFeatureVectors;

julia> x1 = exampleggrv(10, 20);

julia> x2 = exampleggrv(20, 50);

julia> y = vcat(x1, x2);

julia> length(y)
30

julia> length(ranges(y))
70
```
"""
function Base.vcat(A::Vararg{GroupedGenomicRangeVector{T}}) where {T}
    rerange = vcat([ ranges(x) for x in A ]...)
    edata = vcat([ elementdata(x) for x in A ]...)
    recum = copy(A[1].cumulative)

    offset = length(A[1].ranges)
    for i in 2:length(A)
        currange = A[i].ranges
        curcum = A[i].cumulative

        start = length(recum) + 1
        append!(recum, curcum)
        v = view(recum, start:length(recum))
        v .+= offset

        offset += length(A[i].ranges)
    end

    m = combine_metadata(A...)
    return GroupedGenomicRangeVector(rerange, recum, edata, m; islength = false)
end

############# Miscellaneous methods ################

"""
    copy(x::GroupedGenomicRangeVector{T})

Create a copy of `x`.
Note that this function does not copy any of the fields of `x`;
each field in the returned object still refers to the same object as the corresponding field in `x`.

```jldoctest
julia> using GenomicFeatureVectors

julia> x = exampleggrv(5, 10);

julia> y = copy(x);

julia> intervals(y) === intervals(x)
true
```
"""
function Base.copy(x::GroupedGenomicRangeVector{T}) where {T}
    output = GroupedGenomicRangeVector(T)
    output.ranges = x.ranges
    output.cumulative = x.cumulative
    output.elementdata = x.elementdata
    output.metadata = x.metadata
    return output
end

"""
    deepcopy(x::GroupedGenomicRangeVector{T})

Create a deep copy of `x`.

```jldoctest
julia> using GenomicFeatureVectors

julia> x = exampleggrv(5, 10);

julia> y = deepcopy(x);

julia> ranges(y) === ranges(x)
false
```
"""
function Base.deepcopy(x::GroupedGenomicRangeVector{T}) where {T}
    output = GroupedGenomicRangeVector(T)
    output.ranges = deepcopy(x.ranges)
    output.cumulative = deepcopy(x.cumulative)
    output.elementdata = deepcopy(x.elementdata)
    output.metadata = deepcopy(x.metadata)
    return output
end
