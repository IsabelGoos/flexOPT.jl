using DrWatson, JLD2
using SHA
using FilePathsBase: mkpath

# optional parameters to be obtained
_paramget(params::Dict, name::Symbol, default) = get(params, String(name), get(params, name, default))


function hash_parameters(params)
    # Serialize as string and hash
    filtered = Dict(k => v for (k,v) in params)
    str = repr(filtered)
    return bytes2hex(sha1(str))[1:8]  # Short hash
end

"""
    cacheKeyExtras(functionName, paramDict)

Hook for cached functions whose result depends on something that is not
visible in `paramDict` (e.g. the backend actually selected at run time).
Return `nothing` (default: the cache key is exactly the parameters, as before)
or a Dict of extra `String => value` entries; they are added to the parameters
before hashing and stored with the result.
"""
cacheKeyExtras(functionName, paramDict) = nothing

function myProduceOrLoad(functionName,paramDict,directoryName::String)
    output=myProduceOrLoad(functionName,paramDict,directoryName,directoryName)
    return output
end

function myProduceOrLoad(functionName,paramDict,directoryName::String,prefixName::String)

    # this is myFunction strategy for DrWatson

    extras = cacheKeyExtras(functionName, paramDict)
    keyDict = extras === nothing ? paramDict :
        merge(Dict{String,Any}(paramDict), Dict{String,Any}(extras))
    hash_id = hash_parameters(keyDict)


    newDict = Dict{String,Any}(keyDict)
    newDict["hash_id"] = hash_id
    
    output, _ = produce_or_load(functionName,newDict,datadir(directoryName);filename = config -> "$(prefixName)_$(newDict["hash_id"])")
    
    return output

end


# lazy version of myProduceOrLoad

function lazyProduceOrLoad(xString::String, f, args...; folder="./tmp", kwargs...)
    # Ensure the folder exists
    mkpath(folder)
    
    # Build filename

    filename = joinpath(folder, xString * ".jld2")
    
    if isfile(filename)
        # Load saved result
        println("Loading from ", filename)
        @load filename output
        return output
    else
        # Compute and save
        println("Computing ", xString)
        output = isempty(kwargs) ? f(args...) : f(args...; kwargs...)
        @save filename output
        return output
    end
end


function lazyProduceOrLoad(xString::String;folder="./tmp")

    # this is very much lazy
    filename = joinpath(folder, xString * ".jld2")
    if isfile(filename)
        println("This is the laziest ProduceOrLoad")
        println("Loading from ", filename)
        @load filename output
        return output
    else 
        @error(filename, "does not exit! This is the laziest ProduceOrLoad")
        return
    end

end