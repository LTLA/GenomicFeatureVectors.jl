export GenomicRangeVector, intervals, setintervals!
import DataFrames
import GenomicFeatures

############# Class definition ################

"""
The `GenomicRangeVector` is a wrapper around a vector of `GenomicFeatures.Interval` objects,
decorated with some annotations and metadata to make it similar to Bioconductor's `GenomicRanges` objects.
Specifically, we add sequence information (`seqinfo`), per-element data (`elementdata`) and per-vector metadata (`metadata`).
Besides these additions, instances should behave like a vector of intervals.
"""
mutable struct GenomicRangeVector{T} <: GenomicFeatureVector{GenomicFeatures.Interval{T}}
    intervals::Vector{GenomicFeatures.Interval{T}}
    seqinfo::DataFrames.DataFrame
    elementdata::DataFrames.DataFrame
    metadata::Dict{String,Any}

    @doc """
        GenomicRangeVector(T)

    Create an empty `GenomicRangeVector{T}`.

    # Examples
    ```jldoctest
    julia> using GenomicFeatureVectors;

    julia> x = GenomicRangeVector(Nothing);

    julia> typeof(x)
    GenomicRangeVector{Nothing}

    julia> length(x)
    0
    ```
    """
    function GenomicRangeVector(T) 
        mockseq = mock_seqinfo()
        mockelementdata = mock_elementdata(0)
        new{T}(GenomicFeatures.Interval{T}[], mockseq, mockelementdata, Dict{String,Any}())
    end

    @doc """
        GenomicRangeVector(intervals::Vector{GenomicFeatures.Interval{nothing}}))

    Constructor from an existing vector of intervals.

    # Examples
    ```jldoctest
    julia> using GenomicFeatureVectors, GenomicFeatures;

    julia> test = GenomicFeatures.Interval{Nothing}[
               GenomicFeatures.Interval("chrA", 1, 10),
               GenomicFeatures.Interval("chrB", 2, 20),
               GenomicFeatures.Interval("chrC", 3, 30)
           ];

    julia> x = GenomicRangeVector(test);

    julia> length(x)
    3 
    ```
    """
    function GenomicRangeVector(intervals::Vector{GenomicFeatures.Interval{T}}) where {T}
        mockseq = mock_seqinfo()
        mockelementdata = mock_elementdata(length(intervals))
        new{T}(intervals, mockseq, mockelementdata, Dict{String,Any}())
    end

    @doc """
        GenomicRangeVector(
            intervals::Vector{GenomicFeatures.Interval{nothing}}),
            seqinfo::DataFrames.DataFrame, 
            elementdata::DataFrames.DataFrame, 
            metadata::Dict{String,Any} = Dict{String,Any}()
        )

    Full constructor.
    See `setseqinfo!`, `setelementdata!` and `setmetadata!` for details on the requirements of each argument.

    # Examples
    ```jldoctest
    julia> using GenomicFeatureVectors, GenomicFeatures, DataFrames;

    julia> test = GenomicFeatures.Interval{Nothing}[
               GenomicFeatures.Interval("chrA", 1, 10),
               GenomicFeatures.Interval("chrB", 2, 20),
               GenomicFeatures.Interval("chrC", 3, 30)
           ];

    julia> seqinfo = DataFrames.DataFrame(
               seqname = ["chrA", "chrB", "chrC"],
               length = [100, 200, 300],
               circular = [false, false, false],
               genome = ["mm10", "hg19", "rnor6"]
           );

    julia> edata = DataFrames.DataFrame(
               name = ["Alan", "Barney", "Chuck"],
               type = ["protein", "rRNA", "tRNA"]
           );

    julia> meta = Dict("foo" => 2, "bar"=>"BAR");

    julia> x = GenomicRangeVector(test, seqinfo, edata, meta);

    julia> length(x)
    3 
    ```
    """
    function GenomicRangeVector(
            intervals::Vector{GenomicFeatures.Interval{T}}, 
            seqinfo::DataFrames.DataFrame, 
            elementdata::DataFrames.DataFrame, 
            metadata::Dict{String,Any} = Dict{String,Any}()
        ) where {T}

        check_elementdata(elementdata, length(intervals))
        check_seqinfo(seqinfo)
        new{T}(intervals, seqinfo, elementdata, metadata)
    end
end

############# Custom methods ################

"""
    intervals(x::GenomicRangeVector{T})

Retrieve the underlying vector of `GenomicFeatures.Interval` objects from `x`.

# Examples
```jldoctest
julia> using GenomicFeatureVectors

julia> x = examplegrv(10);

julia> y = intervals(x);

julia> length(y)
10

julia> typeof(y[1])
GenomicFeatures.Interval{Nothing}
```
"""
function intervals(x::GenomicRangeVector{T}) where {T}
    return x.intervals
end

function setintervals!(x::GenomicRangeVector{T}, intervals::Vector{GenomicFeatures.Interval{T}}) where {T}
    if length(intervals) != length(x)
        throw(ErrorException("'intervals' and 'x' should have the same length"))
    end
    x.intervals = intervals
    return x
end

############# Vector methods ################

function Base.getindex(x::GenomicRangeVector{T}, I::Int) where {T}
    return x.intervals[I]
end

function Base.getindex(x::GenomicRangeVector{T}, I) where {T}
    return GenomicRangeVector(x.intervals[I], x.seqinfo, x.elementdata[I,:], x.metadata)
end

function Base.setindex!(x::GenomicRangeVector{T}, v::GenomicFeatures.Interval{T}, I::Int) where {T}
    x.intervals[I] = v
    return x
end

function Base.setindex!(x::GenomicRangeVector{T}, v::GenomicRangeVector{T}, I) where {T}
    x.intervals[I] = v.intervals
    x.elementdata[I,:] = v.elementdata
    return x
end

function Base.length(x::GenomicRangeVector{T}) where {T}
    return length(x.intervals)
end

"""
    vcat(A::Vararg{GenomicRangeVector{T}})

Returns a `GenomicRangeVector{T}` containing the concatenated intervals from all `A`.
The length of the output is equal to the sum of the lengths of all input objects.

Element data is concatenated by row, so each entry of `A` should have the same type and name of columns in their [`elementdata`](@ref).

Sequence information in [`seqinfo`](@ref) is concatenated by row across `A`, but only for non-duplicated `seqname` entries.
If duplicates are present, earlier occurrences take precedence and later occurrences are dropped silently.

Metadata in [`metadata`](@ref) is combined across `A`.
Earlier occurrences of each key take precedence, and later occurrences are dropped silently.

# Examples
```jldoctest
julia> using GenomicFeatureVectors

julia> first = examplegrv(5);

julia> second = examplegrv(10);

julia> combined = vcat(first, second);

julia> length(combined)
15
```
"""
function Base.vcat(A::Vararg{GenomicRangeVector{T}}) where {T}
    allint = vcat([intervals(x) for x in A]...)
    alledata = vcat([elementdata(x) for x in A]...)
    allsi = combine_seqinfo(A...)
    allm = combine_metadata(A...)
    return GenomicRangeVector(allint, allsi, alledata, allm)
end

############# Miscellaneous methods ################

"""
    copy(x::GenomicRangeVector{T})

Create a copy of `x`.
Note that this function does not copy any of the fields of `x`;
each field in the returned object still refers to the same object as the corresponding field in `x`.

```jldoctest
julia> using GenomicFeatureVectors

julia> x = examplegrv(10);

julia> y = copy(x);

julia> intervals(y) === intervals(x)
true
```
"""
function Base.copy(x::GenomicRangeVector{T}) where {T}
    output = GenomicRangeVector(T)
    output.intervals = x.intervals
    output.seqinfo = x.seqinfo
    output.elementdata = x.elementdata
    output.metadata = x.metadata
    return output
end

"""
    deepcopy(x::GenomicRangeVector{T})

Create a deep copy of `x`.

```jldoctest
julia> using GenomicFeatureVectors

julia> x = examplegrv(5);

julia> y = deepcopy(x);

julia> intervals(y) === intervals(x)
false
```
"""
function Base.deepcopy(x::GenomicRangeVector{T}) where {T}
    output = GenomicRangeVector(T)
    output.intervals = deepcopy(x.intervals)
    output.seqinfo = deepcopy(x.seqinfo)
    output.elementdata = deepcopy(x.elementdata)
    output.metadata = deepcopy(x.metadata)
    return output
end
