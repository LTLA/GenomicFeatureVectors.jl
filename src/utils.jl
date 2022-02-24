import DataFrames

function check_name(df::DataFrames.DataFrame, msg::String)
    if size(df)[2] < 1 && names(df)[1] != "name"
        throw(ErrorException("'" * string * "' should contain at least one column with the interval names"))
    end

    firstcol = df[!,1]
    if !isa(firstcol, Vector{Nothing}) && !isa(firstcol, AbstractVector{<:AbstractString})
        throw(ErrorException("first column of '" * msg * "' should contain a string vector or nothing"))
    end
end

function check_elementdata(elementdata::DataFrames.DataFrame, expected::Int)
    if expected != size(elementdata)[1]
        throw(DimensionMismatch("number of rows in 'elementdata' is not consistent with vector length (" * string(expected) * ")"))
    end
    check_name(elementdata, "elementadata")
end

function check_seqinfo(seqinfo::DataFrames.DataFrame)
    check_name(seqinfo, "seqinfo")

    if size(seqinfo)[2] < 4
        throw(ErrorException("'seqinfo' should contain at least 4 columns"))
    end

    if names(seqinfo)[2] != "length" || !isa(seqinfo[!,2], AbstractVector{<:Union{Missing,<:Integer}})
        println(names(seqinfo))
        println(seqinfo[!,2])
        throw(ErrorException("second column of 'seqinfo' should be named 'length' and contain integers and/or missing values"))
    end

    if names(seqinfo)[3] != "circular" || !isa(seqinfo[!,3], AbstractVector{<:Union{Missing,Bool}})
        throw(ErrorException("second column of 'seqinfo' should be named 'circular' and contain booleans and/or missing values"))
    end

    if names(seqinfo)[4] != "genome" || !isa(seqinfo[!,4], AbstractVector{<:Union{String,<:AbstractString}})
        throw(ErrorException("second column of 'seqinfo' should be named 'genome' and contain strings and/or missing values"))
    end
end

function mock_seqinfo()
    return DataFrames.DataFrame(
        "name" => String[], 
        "length" => Int[], 
        "circular" => Bool[], 
        "genome" => String[]
    )
end

function mock_elementdata(n::Int)
    return DataFrames.DataFrame("name" => Vector{Nothing}(undef, n))
end

function combine_seqinfo(A::Vararg{GenomicFeatureVector})
    s = copy(seqinfo(A[1]))
    already = Set{String}(s[!,"seqname"])

    for a in 2:length(A)
        curs = seqinfo(A[a])
        curnames = curs[!,"seqname"]

        for i in 1:length(curnames)
            if !(curnames[i] in already)
                push!(already, curnames[i])
                push!(s, curs[i,:])
            end
        end
    end

    return s
end

function combine_metadata(A::Vararg{GenomicFeatureVector})
    m = copy(metadata(A[1]))
    for a in 2:length(A)
        curm = metadata(A[a])
        for (k, v) in curm
            if !haskey(m, k)
                m[k] = v
            end
        end
    end
    return m
end
