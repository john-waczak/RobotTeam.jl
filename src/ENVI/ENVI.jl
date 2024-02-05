# modeled after code found in spectral (SPy)
# https://github.com/spectralpython/spectral/blob/master/spectral/io/envi.py


module ENVI
using HDF5
using DelimitedFiles
using Geodesy
using Dates
using DataInterpolations

export FileNotAnEnviHeader
export EnviHeaderParsingError
export read_envi_header
export get_envi_params
export read_envi_file

export FlightData, HyperspectralImage

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

Base.showerror(io::IO, e::BadENVIHeader) = print(io, "failed to parse ENVI header file.")


# create custom error message for bad header files
struct FileNotAnEnviHeader <: Exception
    file::AbstractString
end

struct EnviHeaderParsingError <: Exception
end

Base.showerror(io::IO, e::FileNotAnEnviHeader) = print(io, e.file, " does not appear to be an ENVI header.")



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

    if inter == "bil" || inter == "BIL"
        img = Array{p["dtype"]}(undef, p["ncols"], p["nbands"], p["nrows"])
        read!(fpath, img)
        img = PermutedDimsArray(img, (2,1,3))
        # img = permutedims(img, (2,1,3))
    elseif inter == "bip" || inter == "BIP"
        img = Array{p["dtype"]}(undef, p["nbands"], p["ncols"], p["nrows"])
        read!(fpath, img)
    else  # BSQ
        img = Array{p["dtype"]}(undef, p["ncols"], p["nrows"], p["nbands"])
        read!(fpath, img)
        img = PermutedDimsArray(img, (3,1,2))
        # img = permutedims(img, (3,1,2))
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











struct FlightData
    start_time
    times
    rolls
    pitches
    headings
    xs
    ys
    zs
    isnorth
    zone
end

nscans(fd::FlightData) = length(fd.times)


"""
    FlightData(lcfpath::String, timespath::String)

Parse `lcfpath` and `timespath` returning a struct containing all relevant flight information for each scanline in the hyper-spectral image.
"""
function FlightData(lcfpath::String, timespath::String; year=2020)
    # 1. compute pixel times
    pixel_ts =  vec(readdlm(timespath, ',', Float64))
    pixel_ts =  pixel_ts .- pixel_ts[1] # assume .lcf and .times start at the same time (t=0.0 s)

    # 2. read flight data
    lcf_data = readdlm(lcfpath, '\t', Float64)

    lcf_ts = @view lcf_data[:,1]
    # start_time = Date.format(gpsToUTC(lcf_ts[1], year), "yyyy-mm-ddTHH:MM:SS.sss")
    # start_time = Date.format(gpsToUTC(lcf_ts[1], year), "yyyy-mm-ddTHH:MM:SS.sss")
    start_time = gpsToUTC(lcf_ts[1], year)
    lcf_times = lcf_ts .- lcf_ts[1]  # so that we start at t=0.0

    roll = -1.0 .*  @view lcf_data[:,2]  # <-- due to opposite convention used by GPS
    #roll_interp = CubicSpline(roll, lcf_times; extrapolate=true)
    roll_interp = LinearInterpolation(roll, lcf_times; extrapolate=true)
    rolls = roll_interp.(pixel_ts)

    pitch = @view lcf_data[:,3]
    #pitch_interp = CubicSpline(pitch, lcf_times; extrapolate=true)
    pitch_interp = LinearInterpolation(pitch, lcf_times; extrapolate=true)
    pitches = pitch_interp.(pixel_ts)

    heading = @view lcf_data[:,4]
    #heading_interp = CubicSpline(heading, lcf_times; extrapolate=true)
    heading_interp = LinearInterpolation(heading, lcf_times; extrapolate=true)
    headings = heading_interp.(pixel_ts)

    lon  = @view lcf_data[:,5]
    lat = @view lcf_data[:,6]
    alt = @view lcf_data[:,7]
    # 3. construct UTMZ coords
    X_lla = LLA.(lat, lon, alt)
    utmzs = [UTMZ(xlla, wgs84) for xlla ∈ X_lla]

    # xs = [utmz.x for utmz ∈ utmzs]
    lcf_x = [utmz.x for utmz ∈ utmzs]
    #x_interp = CubicSpline(lcf_x, lcf_times; extrapolate=true)
    x_interp = LinearInterpolation(lcf_x, lcf_times; extrapolate=true)
    xs = x_interp.(pixel_ts)

    # ys = [utmz.y for utmz ∈ utmzs]
    lcf_y = [utmz.y for utmz ∈ utmzs]
    #y_interp = CubicSpline(lcf_y, lcf_times; extrapolate=true)
    y_interp = LinearInterpolation(lcf_y, lcf_times; extrapolate=true)
    ys = y_interp.(pixel_ts)

    # zs = [utmz.z for utmz ∈ utmzs]
    lcf_z = [utmz.z for utmz ∈ utmzs]
    #z_interp = CubicSpline(lcf_z, lcf_times; extrapolate=true)
    z_interp = LinearInterpolation(lcf_z, lcf_times; extrapolate=true)
    zs = z_interp.(pixel_ts)



    isnorth = utmzs[1].isnorth
    zone = utmzs[1].zone

    return FlightData(start_time, pixel_ts, rolls, pitches, headings, xs, ys, zs, isnorth, zone)
end




struct HyperspectralImage{T0, T1<:AbstractArray}
    # UTMz information
    X::Array{T0, 3}  # x,y,z
    zone::UInt8
    isnorth::Bool

    # Lat/Lon
    Longitudes::Matrix{T0}
    Latitudes::Matrix{T0}

    # Times Information
    Times::Matrix{T0}
    start_time::DateTime

    # Wavelength bins in nm
    λs::Vector{T0}

    # Datacube
    Datacube::T1

    # Viewing Geometry
    Roll::Matrix{T0}
    Pitch::Matrix{T0}
    Heading::Matrix{T0}

    Altitude::Matrix{T0}

    # Camera Geometry
    ViewAngle::Matrix{T0}

    # Solar Geometry
    SolarAzimuth::Matrix{T0}
    SolarElevation::Matrix{T0}
    SolarZenith::Matrix{T0}
end


end
