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

f

hsi = HyperspectralImage(f.bilpath, f.bilhdr, f.lcfpath, f.timespath, f.specpath, f.spechdr; isflipped=true);

hsi_benchmark = @benchmark HyperspectralImage($f.bilpath, $f.bilhdr, $f.lcfpath, $f.timespath, $f.specpath, $f.spechdr; isflipped=true);




# xs, ys, isnorth, zone, Longitudes, Latitudes, IsInbounds, varnames, printnames, λs, Data_μ, Data_σ = resample_datacube(hsi; Δx=Δx)
xs, ys, isnorth, zone, Longitudes, Latitudes, IsInbounds, varnames, printnames, λs, Data_μ = resample_datacube_fast(hsi; Δx=Δx)


resample_benchmark = @benchmark resample_datacube_fast($hsi; Δx=$Δx)

mean(hsi_benchmark)
mean(resample_benchmark)




suite = BenchmarkGroup()

suite["reflectance-conversion"] = BenchmarkGroup()
suite["resampling"] = BenchmarkGroup()

suite["reflectance-conversion"]["Scotty_1-1"] = @benchmarkable HyperspectralImage($f.bilpath, $f.bilhdr, $f.lcfpath, $f.timespath, $f.specpath, $f.spechdr; isflipped=true)

suite["resampling"]["Δx=0.05"] =  @benchmarkable resample_datacube_fast($hsi; Δx=0.05)
suite["resampling"]["Δx=0.1"] =  @benchmarkable resample_datacube_fast($hsi; Δx=0.1)
suite["resampling"]["Δx=0.25"] =  @benchmarkable resample_datacube_fast($hsi; Δx=0.25)
suite["resampling"]["Δx=0.5"] =  @benchmarkable resample_datacube_fast($hsi; Δx=0.5)


tune!(suite);
results = run(suite, verbose=true, seconds=25)


open("./paper/timing-results.txt", "w") do file

    println(file, "---")
    println(file, "Reflectance conversion for Scotty_1-1")
    println(file, "---")
    println(file, mean(results["reflectance-conversion"]["Scotty_1-1"]))
    println(file, "\n")

    println(file, "---")
    println(file, "Resampling results for Δx=0.5")
    println(file, "---")
    println(file, mean(results["resampling"]["Δx=0.5"]))
    println(file, "\n")

    println(file, "---")
    println(file, "Resampling results for Δx=0.25")
    println(file, "---")
    println(file, mean(results["resampling"]["Δx=0.25"]))
    println(file, "\n")

    println(file, "---")
    println(file, "Resampling results for Δx=0.1")
    println(file, "---")
    println(file, mean(results["resampling"]["Δx=0.1"]))
    println(file, "\n")

    println(file, "---")
    println(file, "Resampling results for Δx=0.05")
    println(file, "---")
    println(file, mean(results["resampling"]["Δx=0.05"]))
    println(file, "\n")

end


# create plot of resampling times
sizes = [5, 10, 25, 50]
timings = [
    mean(results["resampling"]["Δx=0.05"]).time / 1e9,
    mean(results["resampling"]["Δx=0.1"]).time / 1e9,
    mean(results["resampling"]["Δx=0.25"]).time / 1e9,
    mean(results["resampling"]["Δx=0.5"]).time / 1e9,
]


my_theme = mints_theme
my_theme.Axis.xticklabelsize=17
my_theme.Axis.yticklabelsize=17
my_theme.Axis.xlabelsize=20
my_theme.Axis.ylabelsize=20
my_theme.Axis.titlesize=22
set_theme!(my_theme)

fig = Figure();
ax = CairoMakie.Axis(fig[1,1], xlabel="grid resolution (cm)", ylabel="execution time (seconds)", title="Hyperspectral Image Reinterpolation")
line  = lines!(ax, sizes, timings, linewidth=3)
scatter  = scatter!(ax, sizes, timings; markersize=15)
fig

save("paper/figures/timing.png", fig)
save("paper/figures/timing.pdf", fig)
save("paper/figures/timing.eps", fig)
save("paper/figures/timing.svg", fig)


