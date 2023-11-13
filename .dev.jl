using RobotTeam
using HDF5
using DelimitedFiles
using Dates

include("utils/config.jl")

# using BenchmarkTools
# using DelimitedFiles
# using Dates

# import CairoMakie as cmk
# using IntervalSets
# using MintsMakieRecipes
# cmk.set_theme!(mints_theme)

# using LoopVectorization
# using LinearAlgebra
# using StaticArrays
# using Meshes, MeshViz

# using Images, ImageContrastAdjustment

# include("utils/vis_tools.jl")



# 1. set up test paths

basepath = "/home/teamlary/gitrepos/utd/RobotTeam.jl/data/raw"
f = get_raw_file_list(get_bil_files(basepath, "NoDye_1")[1])





println("\topening HSI")
hsi = HyperspectralImage(
    f.bilpath,
    f.bilhdr,
    f.lcfpath,
    f.timespath;
    isflipped=true
)





hsi

println("\tresampling to new grid")

xs, ys, isnorth, zone, Longitudes, Latitudes, IsInbounds, varnames, printnames, λs, Data = resample_datacube(hsi)

println("\tconverting to reflectance")
generateReflectance!(Data, f.specpath, f.spechdr, λs)

println("\tgenerating derived metrics")
generate_derived_metrics!(Data, IsInbounds, varnames, λs)

fname = split(f.lcfpath, "/")[end-1] * ".h5"

println("\tsaving as $fname")
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
    joinpath(outpath, day, run, fname);
    Δx=Δx,
    is_spec_chunked=is_spec_chunked,
    is_band_chunked=is_band_chunked
)



