function check_name(df::DataFrames.DataFrame, msg::String)
    if size(mcols)[2] < 1 && names(mcols)[1] != "name"
        throw(ErrorException("'" * string "' should contain at least one column with the interval names"))
    end

    firstcol = mcols[!,1]
    if !isa(firstcol, Vector{Nothing}) && !isa(firstcol, AbstractVector{<:AbstractString})
        throw(ErrorException("first column of '" * msg "' should contain a string vector or nothing"));
    end
end

function check_mcols(mcols::DataFrames.DataFrame, expected::Int)
    if expected != size(mcols)[1]
        throw(DimensionMismatch("number of rows in 'mcols' is not consistent with vector length (" * string(expected) * ")"))
    end
    check_name(mcols, "mcols")
end

function check_seqinfo(seqinfo::DataFrames.DataFrame)
    check_name(mcols, "seqinfo")

    if size(seqinfo)[2] < 4
        throw(ErrorException("'seqinfo' should contain at least 4 columns"))
    end

    if names(seqinfo)[2] != "length" || !isa(seqinfo[!,2], AbstractVector{Integer})
        throw(ErrorException("second column of 'seqinfo' should be named 'lengths' and contain integer lengths"))
    end

    if names(seqinfo)[3] != "circular" || !isa(seqinfo[!,2], AbstractVector{Bool})
        throw(ErrorException("second column of 'seqinfo' should be named 'circular' and contain circular flags"))
    end

    if names(seqinfo)[4] != "genome" || !isa(seqinfo[!,2], AbstractVector{AbstractString})
        throw(ErrorException("second column of 'seqinfo' should be named 'genome' and contain genome identifiers"))
    end
end

function mock_seqinfo()
    return DataFrames.DataFrame(name = String[], length=Int[], circular=Bool[], genome=String[])
end

function mock_mcols(n::Int)
    return DataFrames.DataFrame(name = Vector{Nothing}(undef, n))
end


