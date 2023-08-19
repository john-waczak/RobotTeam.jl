using RobotTeam
using HDF5
using BenchmarkTools
import CairoMakie as cmk
using MintsMakieRecipes
cmk.set_theme!(mints_theme)

include("utils/vis_tools.jl")
include("utils/config.jl")

basepath

files = get_raw_file_list.(get_bil_files(basepath, "Scotty_1"))

f = files[1]

hsi = HyperspectralImage(f.bilpath, f.bilhdr, f.lcfpath, f.timespath, f.specpath, f.spechdr; isflipped=true);


xs, ys, isnorth, zone, Longitudes, Latitudes, IsInbounds, varnames, printnames, λs, Data_μ, Data_σ = resample_datacube(hsi; Δx=Δx)


idxs = findall(IsInbounds)
Data_μ[1,idxs[1]]
eltype(Data_μ)




t = @benchmark HyperspectralImage(f.bilpath, f.bilhdr, f.lcfpath, f.timespath, f.specpath, f.spechdr; isflipped=true);
tune!(t)


# make sure to run the timings for multiple values of Δx e.g.  0.05, 0.1, 0.25, 0.5

# fib(n) = n <= 1 ?  1 : fib(n - 2) + fib(n - 1)
# suite = BenchmarkGroup()
# suite["fib"] = BenchmarkGroup(["tag1", "tag2"])
# suite["fib"][10] = @benchmarkable fib(10)
# suite["fib"][20] = @benchmarkable fib(20)

# t = @benchmark fib(10)

# tune!(suite)
# results = run(suite, verbose = true)

# BenchmarkTools.save("output.json", median(results))



