using RobotTeam
using BenchmarkTools
using HDF5
using CairoMakie
using MintsMakieRecipes

set_theme!(mints_theme)

basepath = "/home/teamlary/gitrepos/utd/RobotTeam.jl/data/raw"

# three HSI files to test
f1 = get_raw_file_list(get_bil_files(basepath, "Dye_1")[1])
f2 = get_raw_file_list(get_bil_files(basepath, "NoDye_1")[1])
f3 = get_raw_file_list(get_bil_files(basepath, "Scotty_1")[1])



@benchmark read_envi_header($f1.bilhdr)  # 107 μs

h = read_envi_header(f1.bilhdr)

@benchmark get_envi_params($h)  # 46 μs

params = get_envi_params(h)


# this reads in the relevant flight data time-stamped to the correct
@benchmark FlightData($f1.lcfpath, $f1.timespath)  # 3.126 ms
fd = FlightData(f1.lcfpath, f1.timespath)
fd.start_time



@benchmark read_envi_file($f1.bilpath, $f1.bilhdr)  # 337 ms
@benchmark read_envi_file($f1.specpath, $f1.spechdr)  # 651.675 μs





# 1. read the raw radiance cube

Datacube, h, p = read_envi_file(f1.bilpath, f1.bilhdr)

size(Datacube)

X = zeros(3, size(Datacube,2), size(Datacube,3))
Latitudes  = zeros(size(Datacube,2), size(Datacube,3))
Longitudes = zeros(size(Datacube,2), size(Datacube,3))
Roll = zeros(size(Datacube,2), size(Datacube,3))
Pitch = zeros(size(Datacube,2), size(Datacube,3))
Heading = zeros(size(Datacube,2), size(Datacube,3))
Times = zeros(size(Datacube,2), size(Datacube,3))
ViewAngle = zeros(size(Datacube,2), size(Datacube,3))
SolarAzimuth = zeros(size(Datacube,2), size(Datacube,3))
SolarElevation = zeros(size(Datacube,2), size(Datacube,3))
SolarZenith = zeros(size(Datacube,2), size(Datacube,3))
fd.start_time

@benchmark generateCoords!(
    X,
    Longitudes,
    Latitudes,
    Roll,
    Pitch,
    Heading,
    Times,
    ViewAngle,
    SolarAzimuth,
    SolarElevation,
    SolarZenith,
    fd.start_time,
    fd
)  # 692 ms !!!


generateCoords!(
    X,
    Longitudes,
    Latitudes,
    Roll,
    Pitch,
    Heading,
    Times,
    ViewAngle,
    SolarAzimuth,
    SolarElevation,
    SolarZenith,
    fd.start_time,
    fd
)  # 45.283 ms w/ no lat/lon conversion 1.045 s

typeof(Datacube) <: Array

hsi = HyperspectralImage(f1.bilpath, f1.bilhdr, f1.lcfpath, f1.timespath; isflipped=true)
hsi = HyperspectralImage(f1.bilpath, f1.bilhdr, f1.lcfpath, f1.timespath)


@benchmark HyperspectralImage($f1.bilpath, $f1.bilhdr, $f1.lcfpath, $f1.timespath)


@benchmark resample_datacube($hsi)  # 3.5 seconds


xs, ys, isnorth, zone, Longitudes, Latitudes, IsInbounds, varnames, printnames, λs, Data = resample_datacube(hsi)


@benchmark generateReflectance!(Data, f1.specpath, f1.spechdr, λs) # 104 ms
@benchmark generate_derived_metrics!(Data, IsInbounds, varnames, λs)  # 250 ms




generateReflectance!(Data, f1.specpath, f1.spechdr, λs) # 104 ms
generate_derived_metrics!(Data, IsInbounds, varnames, λs)  # 250 ms


save_resampled_hsi(
    xs,
    ys,
    isnorth,
    zone,
    Longitudes,
    Latitudes,
    IsInbounds,
    varnames,
    printnames,
    λs,
    Data,
    "./test.h5"
)


ij_pixels = findall(IsInbounds)

rgb = getRGB(Data, λs, ij_pixels)


heatmap(Longitudes, Latitudes, rgb)

# verify that the color looks right
