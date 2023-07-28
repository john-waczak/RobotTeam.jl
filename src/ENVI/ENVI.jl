# modeled after code found in spectral (SPy)
# https://github.com/spectralpython/spectral/blob/master/spectral/io/envi.py


module ENVI
using HDF5
using DelimitedFiles

export FileNotAnEnviHeader
export EnviHeaderParsingError
export read_envi_header
export get_envi_params
export read_envi_file
export envi_to_hdf5


envi_to_dtype = Dict(
    "1" => UInt8,
    "2" => Int16,
    "3" => Int32,
    "4" => Float32,
    "5" => Float64,
    "6" => ComplexF32,
    "9" => ComplexF64,
    "12" => UInt16,
    "13" => UInt32,
    "14" => Int64,
    "15" => UInt64,
)

dtype_to_envi = Dict(val => key for (key,val) ∈ envi_to_dtype)



struct BadENVIHeader <: Exception
    fname::AbstractString
end

Base.showerror(io::IO, e::BadENVIHeader) = print(io, e.fname, " does not appear to be an ENVI header.")
Base.showerror(io::IO, e::BadENVIHeader) = print(io, "failed to parse ENVI header file.")



# create custom error message for bad header files
struct FileNotAnEnviHeader <: Exception
    file::AbstractString
end

struct EnviHeaderParsingError <: Exception
end

Base.showerror(io::IO, e::FileNotAnEnviHeader) = print(io, e.file, " does not appear to be an ENVI header.")
Base.showerror(io::IO, e::FileNotAnEnviHeader) = print(io, "failed to parse ENVI header file.")



"""
    read_envi_header(file::String)

Reads and ENVI `.hdr` file header and returns the parameters in a dictionary as strings.
"""
function read_envi_header(file::String)

    f = open(file, "r")
    # make sure we have a header file by checking first line
    starts_with_ENVI = startswith(readline(f), "ENVI")
    if !(starts_with_ENVI)
        throw(FileNotAnEnviHeader(file))
    end

    lines = readlines(f)
    close(f)

    res = Dict()

    try
        for line ∈ lines
            if occursin("=", line) && line[1] != ';'
                splitline = split(line, "=")
                key = strip(splitline[1])
                val = strip(splitline[2])

                res[key] = val

                # check for array information
                if val[1] == '{'
                    if key == "description"
                        res[key] = strip(val[2:end-1])
                    else
                        vals = [strip(v) for v ∈ split(val[2:end-1], ",")]
                        res[key] = vals
                    end
                else
                    res[key] = val
                end

            end
        end
        return res
    catch e
        throw(EnviHeaderParsingError)
    end


    # make sure we have mandatory parameters
    mandatory_params = ["lines", "samples", "bands", "data type", "interleave", "byte order"]
    if any([!(mp ∈ keys(res)) for mp ∈ mandatory_params])
        throw(EnviHeaderParsingError, "Missing at least one mandatory parameter")
    end

    return res
end




"""
    get_envi_params(h::Dict)

Parse dict returned by `read_envi_Header` and return parameters needed for reading binary file.
"""
function get_envi_params(h::Dict)
    params = Dict()

    params["header offset"] = 0
    for (key, val) ∈ h
        if key == "wavelength"
            params["wavelengths"] = parse.(Float64, val)
        elseif key == "timestamp"
            params[key] = val
        elseif key == "interleave"
            params[key] = val
        elseif key == "shutter"
            params[key] = parse(Float64, val)
        elseif key == "header offset"
            params[key] = parse(Int, val)

        elseif key == "bands"
            params["nbands"] = parse(Int, val)
        elseif key == "lines"
            params["nrows"] = parse(Int, val)
        elseif key == "samples"
            params["ncols"] = parse(Int, val)
        elseif key == "header offset"
            params["offset"] = parse(Int, val)
        elseif key == "byte order"
            params[key] = parse(Int, val)
        elseif key == "data type"
            params["dtype"] = envi_to_dtype[val]
        else
            params[key] = val
        end
    end

    return params
end



"""
    read_envi_file(fpath::String, hdrpath::String)

Read an ENVI formatted HSI file located at `fpath` with it's associated metatdata in `hdrpath`. Returns an image array `img`, the parsed header dictionary and the parameter dictionary.
"""
function read_envi_file(fpath::String, hdrpath::String)
    # read header file
    h = read_envi_header(hdrpath)
    p = get_envi_params(h)
    inter = h["interleave"]

    # change to bip for all as it will be most memory efficient when looping
    if inter == "bil" || inter == "BIL"
        img = Array{p["dtype"]}(undef, p["ncols"], p["nbands"], p["nrows"])
        read!(fpath, img)
        # img = PermutedDimsArray(img, (2,1,3))
        img = permutedims(img, (2,1,3))
    elseif inter == "bip" || inter == "BIP"
        img = Array{p["dtype"]}(undef, p["nbands"], p["ncols"], p["nrows"])
        read!(fpath, img)
    else
        img = Array{p["dtype"]}(undef, p["ncols"], p["nrows"], p["nbands"])
        # img = PermutedDimsArray(img, (3,1,2))
        img = permutedims(img, (3,1,2))
    end

    h["interleave"] = "bip"
    h["shape"] = "(band,col,row)"
    return img, h, p
end




"""
    envi_to_hdf5(fpath::String, hdrpath::String, lcfpath::outpath::String)

Read ENVI formatted HSI file from `fpath` and its associated metadata file `hdrpath`. Save the array to an hdf5 file at `outapth`.
"""
function envi_to_hdf5(
    bilpath::String,
    bilhdr::String,
    lcfpath::String,
    timespath::String,
    specpath::String,
    spechdr::String,
    outpath::String,
    )


    h5open(outpath, "cw") do fid
        g = create_group(fid, "raw")

        # create subgroups for each data set type
        rad = create_group(g, "radiance")
        down = create_group(g, "downwelling")
        lcf = create_group(g, "lcf")
        times = create_group(g, "times")

        # use let block to keep img from persisting
        let
            img, h, p = read_envi_file(bilpath, bilhdr)

            # write the envi radiance data
            rad["radiance", chunk=(p["nbands"], 1, 1)] = img
            # rad["wavelengths"] = parse.(Float64, h["wavelength"])
            for (key, val) ∈ p
                if key != "dtype"
                    rad[key] = val
                end
            end
        end

        let
            spec, hspec, pspec = read_envi_file(specpath, spechdr)
            # write the downwelling irradiance
            down["irradiance"] = spec
            for (key, val) ∈ pspec
                if key != "dtype"
                    down[key] = val
                end
            end
        end

        lcf["lcf"] = readdlm(lcfpath, '\t', Float64)

        ts =  readdlm(timespath, ',', Float64)

        times["times"] = ts .- ts[1]  # assume .lcf and .times start at the same time
    end
end


end
