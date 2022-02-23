export GenomicRangeVector
import DataFrames
import GenomicFeatures

struct GenomicRangeVector <: GenomicFeatureVector
    intervals::Vector{GenomicFeatures.Interval{nothing}}
    seqinfo::DataFrames.DataFrame
    mcols::DataFrames.DataFrame
    metadata::Dict{String,Any}

    function GenomicRangeVector()
        mockseq = mock_seqinfo()
        mockmcols = mock_mcol(length(intervals))
        new(GenomicFeatures.Interval{nothing}[], mockseq, mockmcols, Dict{String,Any}())
    end

    function GenomicRangeVector(intervals::Vector{GenomicFeatures.Interval{nothing}})
        mockseq = mock_seqinfo()
        mockmcols = mock_mcol(length(intervals))
        new(intervals, mockseq, mockmcols, Dict{String,Any}())
    end

    function GenomicRangeVector(
            intervals::Vector{GenomicFeatures.Interval{nothing}}, 
            seqinfo::DataFrames.DataFrame, 
            mcols::DataFrames.DataFrame, 
            metadata::Dict{String,Any} = Dict{String,Any}()
        )

        check_mcols(mcols, length(intervals))
        check_seqinfo(seqinfo)
        new(intervals, seqinfo, mcols, metadata)
    end
end

function intervals(x::GenomicRangeVector)
    return x.intervals
end

function setintervals!(x::GenomicRangeVector, intervals::Vector{GenomicFeatures.Interval{nothing}})
    if length(intervals) != length(x)
        throw(ErrorException("'intervals' and 'x' should have the same length"))
    end
    x.intervals = intervals
    return x
end

function mcols(x::GenomicRangeVector)
    return x.mcols
end

function setmcols!(x::GenomicRangeVector, mcols::DataFrames.DataFrame)
    if size(mcols)[1] != length(x)
        throw(ErrorException("'intervals' and 'x' should have the same length"))
    end
    x.mcols = mcols
    return x
end

function metadata(x::GenomicRangeVector)
    return x.metadata
end

function setmetadata!(x::GenomicRangeVector, metadata::Dict{String,Any})
    x.metadata = metadata
    return x
end

function Base.getindex(x::GenomicRangeVector, I::Int)
    return x.intervals[I]
end

function Base.getindex(x::GenomicRangeVector, I)
    return GenomicRangeVector(x.intervals[I], x.mcols[I,:], x.metadata)
end

function Base.setindex!(x::GenomicRangeVector, v::GenomicFeatures.Interval{nothing}, I::Int)
    x.intervals[I] = v
    return x
end

function Base.setindex!(x::GenomicRangeVector, v::GenomicRangeVector, I)
    x.intervals[I] = v.intervals
    x.mcols[I,:] = v.mcols
    x.metadata = v.metadata
    return x
end

function Base.length(x::GenomicRangeVector)
    return length(x.intervals)
end
