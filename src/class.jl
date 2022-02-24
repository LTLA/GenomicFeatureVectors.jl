export GenomicFeatureVector
export elementdata, seqinfo, metadata
export setelementdata!, setseqinfo!, setmetadata!
import DataFrames

abstract type GenomicFeatureVector{T} <: AbstractVector{T} end

function Base.iterate(x::GenomicFeatureVector{T}) where {T}
    if length(x) > 0
        return (x[1], 1)
    else
        return nothing
    end
end

function Base.iterate(x::GenomicFeatureVector{T}, state::Int) where {T}
    if state >= length(x)
        return nothing
    else
        state += 1
        return (x[state], state)
    end
end

function Base.size(x::GenomicFeatureVector{T}) where {T}
    return (length(x),)
end

function Base.size(x::GenomicFeatureVector{T}, dim::Int) where {T}
    # Assuming 'dim' is 1.
    return length(x)
end

function seqinfo(x::GenomicFeatureVector{T}) where {T}
    return x.seqinfo
end

function setseqinfo!(x::GenomicFeatureVector{T}, seqinfo::DataFrames.DataFrame) where {T}
    check_seqinfo(seqinfo)
    return x.seqinfo
end

function elementdata(x::GenomicFeatureVector{T}) where {T}
    return x.elementdata
end

function setelementdata!(x::GenomicFeatureVector{T}, elementdata::DataFrames.DataFrame) where {T}
    check_elementdata(elementdata, length(x))
    x.elementdata = elementdata
    return x
end

function metadata(x::GenomicFeatureVector{T}) where {T}
    return x.metadata
end

function setmetadata!(x::GenomicFeatureVector{T}, metadata::Dict{String,Any}) where {T}
    x.metadata = metadata
    return x
end
