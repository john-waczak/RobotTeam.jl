# modeled after code found in spectral (SPy)
# https://github.com/spectralpython/spectral/blob/master/spectral/io/envi.py


module ENVI
using HDF5
using DelimitedFiles
using Geodesy
using Dates


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
    function leap_count(year::Int)

Determine the number of leap seconds introduced before the given date. See the links below:
- [Stack Overflow Post](https://stackoverflow.com/questions/33415475/how-to-get-current-date-and-time-from-gps-unsegment-time-in-python)
- [Official Leap Seconds Table](https://hpiers.obspm.fr/eop-pc/index.php?index=TAI-UTC_tab&lang=en)
"""
function leap_count(year::Int)
    leap_years =[1980,
                 1981,
                 1982,
                 1983,
                 1985,
                 1988,
                 1990,
                 1991,
                 1992,
                 1993,
                 1994,
                 1996,
                 1997,
                 1999,
                 2006,
                 2009,
                 2012,
                 2015,
                 2017,
                 ];

    leap_secs = [19,
                 20,
                 21,
                 22,
                 23,
                 24,
                 25,
                 26,
                 27,
                 28,
                 29,
                 30,
                 31,
                 32,
                 33,
                 34,
                 35,
                 36,
                 37,];


    idx = findfirst(year .< leap_years)
    if idx == nothing
        return leap_secs[end]
    else
        return leap_secs[idx - 1]
    end

end


"""
    function gpsToUTC(gps_sec, year)

Given the gps_time in seconds and the current year, return the UTC time accounting for leap seconds. See this [Stack Overflow Post](https://stackoverflow.com/questions/33415475/how-to-get-current-date-and-time-from-gps-unsegment-time-in-python) for more details.
"""
function gpsToUTC(gps_sec, year)
    sec = round(gps_sec)
    ms = (gps_sec - sec)
    # convert to DateTime object
    sec = Second(Int(sec))
    ms = Millisecond(Int(round(1000*ms, digits=0)))

    gps = sec + ms

    utc = DateTime(1980, 1, 6) + (gps - Second(leap_count(year) - leap_count(1980)) )
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
    outpath::String;
    )

    # assume year is 2019 and overwrite with year from hdr file
    year = 2019

    h5open(outpath, "cw") do fid
        println("\tcreating groups")
        g = create_group(fid, "raw")

        # create subgroups for each data set type
        rad = create_group(g, "radiance")
        down = create_group(g, "downwelling")
        lcf = create_group(g, "lcf")
        times = create_group(g, "times")

        # use let block to keep img from persisting
        let
            println("\treading radiance")
            img, h, p = read_envi_file(bilpath, bilhdr)

            year = parse(Int, split(split(p["timestamp"], "-")[1], "/")[end])
            println("\tyear taken: $(year)")

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
            println("\treading downwelling irradiance")
            spec, hspec, pspec = read_envi_file(specpath, spechdr)
            # write the downwelling irradiance
            down["irradiance"] = spec
            for (key, val) ∈ pspec
                if key != "dtype"
                    down[key] = val
                end
            end
        end

        println("\treading flight data")
        lcf_data = readdlm(lcfpath, '\t', Float64)

        lon  = lcf_data[:,5]
        lat = lcf_data[:,6]
        alt = lcf_data[:,7]

        lcf["longitude"] = lon
        lcf["latitude"] = lat
        lcf["altitude"] = alt

        # generate x,y,z positions in
        println("\tconstructing local coordinates")
        X_lla = LLA.(lat, lon, alt)
        utmzs = [UTMZ(xlla, wgs84) for xlla ∈ X_lla]

        # position in meters in wgs84 UTMZ ellipsoid
        lcf["x"] = [utmz.x for utmz ∈ utmzs]
        lcf["y"] = [utmz.y for utmz ∈ utmzs]
        lcf["z"] = [utmz.z for utmz ∈ utmzs]
        lcf["isnorth"] = [utmz.isnorth for utmz ∈ utmzs]
        lcf["zone"] = [utmz.zone for utmz ∈ utmzs]


        # apply corrections to orientation data
        utm_zones = range(-180, stop=180, step=6)  # utm zones are every 6 degrees

        # heading adjustment due to convergence of lines of longitude
        # towards the poles
        zones = [utm_zones[utmz.zone] for utmz ∈ utmzs]
        Δ = atan.(tand.(lon .- (zones .+ 3.0)).*sind.(lat))

        heading_correct = lcf_data[:,4]  .- Δ  # <-- this is what's used in the paper by Muller
        # final assignments
        lcf["roll"] = -lcf_data[:,2]  # <-- due to opposite convention used by GPS
        lcf["pitch"] = lcf_data[:,3]
        lcf["heading"] = heading_correct


        # now we compute the times as measured by the IMU in seconds
        # NOTE: this sensor uses the weird gps time which measures
        # seconds since Jan 6th 1980 without accounting for leap year
        # see above function notes for more details.
        ts = lcf_data[:,1]
        lcf["start-time"] = Dates.format(gpsToUTC(ts[1], year), "yyyy-mm-ddTHH:MM:SS.sss")
        lcf["times"] = ts .- ts[1]  # so that we start at t=0.0


        # now save corresponding scanline times from the .times file
        # these come from the camera, not the IMU
        println("\tsaving times")
        ts =  readdlm(timespath, ',', Float64)

        times["times"] = ts .- ts[1]  # assume .lcf and .times start at the same time
    end
end




end
