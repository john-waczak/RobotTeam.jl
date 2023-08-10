module Georectification

using ..ENVI: read_envi_header, get_envi_params, FlightData, HyperspectralImage

using HDF5
using StaticArrays
using LinearAlgebra
using Geodesy
using Meshes
using RelocatableFolders
using SolarGeometry

# set up calibration path
const calibration_path = @path normpath(joinpath(@__DIR__, "../../assets", "calibration"))
const gain_path = @path normpath(joinpath(@__DIR__, "../../assets", "calibration", "gain.spec"))
const gain_hdr = @path normpath(joinpath(@__DIR__, "../../assets", "calibration", "gain.spec.hdr"))
const offset_path = @path normpath(joinpath(@__DIR__, "../../assets", "calibration", "offset.spec"))
const offset_hdr = @path normpath(joinpath(@__DIR__, "../../assets", "calibration", "offset.spec.hdr"))

@assert ispath(calibration_path)
@assert ispath(gain_path)
@assert ispath(gain_hdr)
@assert ispath(offset_path)
@assert ispath(offset_hdr)

include("reflectance.jl")
export generateReflectance!
export generateDerivedMetrics!

include("georectify.jl")
export generateCoords!


include("mesh.jl")
export build_mesh


include("resample.jl")
export bump_to_nearest_Δx, get_new_bounds, get_resampled_grid



function HyperspectralImage(
    bilpath::String,
    bilhdr::String,
    lcfpath::String,
    timespath::String,
    specpath::String,
    spechdr::String;
    θ_view=30.8,
    z_ground=292.0,
    isflipped=false
    )
    # 0. load flight data
    fdata = FlightData(lcfpath, timespath)  # 5.861 ms

    # 1. read header
    h = read_envi_header(bilhdr)
    p = get_envi_params(h)

    nbands = p["nbands"]
    nsamples = p["ncols"]
    nscanlines = p["nrows"]

    # 2. preallocate outputs
    X = Array{Float64}(undef, 3, nsamples, nscanlines);
    zone = fdata.zone
    isnorth = fdata.isnorth
    Longitudes = Matrix{Float64}(undef, nsamples, nscanlines);
    Latitudes = Matrix{Float64}(undef, nsamples, nscanlines);
    Times = Matrix{Float64}(undef, nsamples, nscanlines);
    start_time = fdata.start_time
    λs = p["wavelengths"]

    Reflectance = Array{Float64}(undef, nbands, nsamples, nscanlines)

    Roll = Matrix{Float64}(undef, nsamples, nscanlines);
    Pitch = Matrix{Float64}(undef, nsamples, nscanlines);
    Heading = Matrix{Float64}(undef, nsamples, nscanlines);

    ViewingAngle = Matrix{Float64}(undef, nsamples, nscanlines);

    SolarAzimuth = Matrix{Float64}(undef, nsamples, nscanlines);
    SolarElevation = Matrix{Float64}(undef, nsamples, nscanlines);
    SolarZenith = Matrix{Float64}(undef, nsamples, nscanlines);


    # 3. Compute Reflectance Data
    let
        Radiance, h, p = read_envi_file(bilpath, bilhdr)
        generateReflectance!(Reflectance, Radiance, specpath, spechdr, λs)
    end

    # 3. Generate Coordinates
    generateCoords!(X,Longitudes, Latitudes, Roll,Pitch,Heading,Times,ViewingAngle,SolarAzimuth,SolarElevation,SolarZenith,start_time,fdata;θ_view=θ_view,z_ground=z_ground,isflipped=isflipped)


    return HyperspectralImage(
        X,
        zone,
        isnorth,
        Longitudes,
        Latitudes,
        Times,
        start_time,
        λs,
        Reflectance,
        Roll,
        Pitch,
        Heading,
        ViewingAngle,
        SolarAzimuth,
        SolarElevation,
        SolarZenith
    )
end





end
