export MultiGenomicRangeVector
export intervals, setintervals!

function create_cumulative(groups::AbstractVector{Integer}, islength::Bool)
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
                if i > 1 && cumulative[i] < cumulative[i-1]
                    throw(ErrorException("endpoints in 'groups' must be non-decreasing"))
                end
                cumulative[i] = groups[i]
            end
        end
    end

    return cumulative
end

struct MultiGenomicRangeVector{T} <: GenomicFeatureVector{GenomicRangeVector{T}}
    ranges::GenomicRangeVector{T}
    cumulative::Vector{Int}
    elementadata::DataFrames.DataFrame
    metadata::Dict{String,Any}

    function MultiGenomicRangeVector() where {T}
        mockelementdata = mock_elementdata(length(intervals))
        new{T}(GenomicRangeVector{T}(), Int[], mockelementdata, Dict{String,Any})
    end

    function MultiGenomicRangeVector(ranges::GenomicRangeVector{T}, groups::Vector{Int}; islength = true) where {T}
        cumulative = create_cumulative(groups, islength)
        mockelementdata = mock_elementdata(length(intervals))
        new{T}(ranges, cumulative, mockelementdata, Dict{String,Any})
    end

    function MultiGenomicRangeVector(
            ranges::GenomicRangeVector{T}, 
            groups::Vector{Int}, 
            elementadata::DataFrames.DataFrame, 
            metadata::Dict{String,Any}; 
            islength = true
        ) where {T}

        cumulative = create_cumulative(groups, islength)
        new{T}(ranges, cumulative, elementdata, metadata)
    end
end

function length(x::MultiGenomicRangeVector{T}) where {T}
    return length(x.cumulative)
end

function intervals(x::MultiGenomicRangeVector{T}) where {T}
    return intervals(x.ranges)
end

function setintervals!(x::MultiGenomicRangeVector{T}, intervals::GenomicRangeVector{T}) where {T}
    setintervals!(x, intervals)
    return 
end

function ranges(x::MultiGenomicRangeVector{T}) where {T}
    return x.ranges
end

function setranges!(x::MultiGenomicRangeVector{T}, value::GenomicRangeVector{T}) where {T}
    if length(value) != length(x.ranges)
        throw(ErrorException("lengths of 'value' and 'ranges(x)' must be the same"))
    end
    x.ranges = value
    return x
end

function lengths(x::MultiGenomicRangeVector{T}) where {T}
    output = Vector{Int}(undef, length(x))
    if length(output) > 0
        output[1] = x.cumulative[1]
        for i in 2:length(output)
            output[i] = x.cumulative[i] - start[i-1]
        end
    end
    return output
end

function single_element_range(x::MultiGenomicRangeVector{T}, i::Int) where {T}
    s = (i == 1 ? 0 : cumulative[1]) + 1
    l = start[i]
    return s:l
end

function Base.getindex(x::MultiGenomicRangeVector{T}, i::Int) where {T}
    rr = single_element_range(x, i)
    return x.ranges[rr]
end

function Base.setindex!(x::MultiGenomicRangeVector{T}, value::GenomicRangeVector{T}, i::Int) where {T}
    rr = single_element_range(x, i)
    if length(value) == length(rr)
        x.ranges[rr] = value
    else
        # Reassembly required if lengths don't match up.
        before = 1:(first(rr) - 1)
        after = (last(rr) + 1):length(x.ranges)
        x.ranges = vcat(x.ranges[before], value, x.ranges[after])
        delta = length(rr) - length(value)
        for j in i:length(x.cumulative)
            x.cumulative[j] += delta
        end
    end
    return x
end

function Base.getindex(x::MultiGenomicRangeVector{T}, I) where {T}
    idx = to_indices(I, x)

    new_range_index = Int[]
    latest = 0
    new_cum = Int[]
    for i in idx
        rr = single_element_range(x, i)
        append!(new_range_index, rr)
        latest += length(rr)
        push!(new_cum, latest)
    end

    return MultiGenomicRangeVector(x.ranges[new_range_index], new_cum, x.elementdata[idx,:], x.metadata; islength = false)
end
