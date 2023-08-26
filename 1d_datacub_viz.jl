using RobotTeam
using HDF5
using BenchmarkTools
using CairoMakie
using MintsMakieRecipes

set_theme!(mints_theme)

include("utils/vis_tools.jl")
include("utils/config.jl")

basepath

files = get_raw_file_list.(get_bil_files(basepath, "Scotty_1"))

f = files[1]


hsi = HyperspectralImage(f.bilpath, f.bilhdr, f.lcfpath, f.timespath, f.specpath, f.spechdr; isflipped=true);
xs, ys, isnorth, zone, Longitudes, Latitudes, IsInbounds, varnames, printnames, λs, Data_μ = resample_datacube_fast(hsi; Δx=Δx)



# plot pre-georectified HSI

# plot georectified HSI

# plot HSI on map

