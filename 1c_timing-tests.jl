using RobotTeam
using HDF5
using BenchmarkTools
using CairoMakie
using MintsMakieRecipes

set_theme!(mints_theme)

include("utils/vis_tools.jl")
include("utils/config.jl")


basepath = "/home/teamlary/gitrepos/utd/RobotTeam.jl/data/raw"
f1 = get_raw_file_list(get_bil_files(basepath, "Dye_1")[1])
f2 = get_raw_file_list(get_bil_files(basepath, "NoDye_1")[1])
f3 = get_raw_file_list(get_bil_files(basepath, "Scotty_1")[1])

f = f3


hsi = HyperspectralImage(f.bilpath, f.bilhdr, f.lcfpath, f.timespath, f.specpath, f.spechdr; isflipped=true);


suite = BenchmarkGroup()

suite["reflectance-conversion"] = BenchmarkGroup()
suite["resampling"] = BenchmarkGroup()

suite["reflectance-conversion"]["Dye_1"] = @benchmarkable HyperspectralImage($f1.bilpath, $f1.bilhdr, $f1.lcfpath, $f1.timespath, $f1.specpath, $f1.spechdr; isflipped=true)
suite["reflectance-conversion"]["NoDye_1"] = @benchmarkable HyperspectralImage($f2.bilpath, $f2.bilhdr, $f2.lcfpath, $f2.timespath, $f2.specpath, $f2.spechdr; isflipped=true)
suite["reflectance-conversion"]["Scotty_1"] = @benchmarkable HyperspectralImage($f3.bilpath, $f3.bilhdr, $f3.lcfpath, $f3.timespath, $f3.specpath, $f3.spechdr; isflipped=true)


suite["resampling"]["Δx=0.05"] = @benchmarkable resample_datacube_fast($hsi; Δx=0.05)
suite["resampling"]["Δx=0.1"] = @benchmarkable resample_datacube_fast($hsi; Δx=0.1)
suite["resampling"]["Δx=0.2"] = @benchmarkable resample_datacube_fast($hsi; Δx=0.2)
suite["resampling"]["Δx=0.3"] = @benchmarkable resample_datacube_fast($hsi; Δx=0.3)
suite["resampling"]["Δx=0.4"] = @benchmarkable resample_datacube_fast($hsi; Δx=0.4)
suite["resampling"]["Δx=0.5"] =  @benchmarkable resample_datacube_fast($hsi; Δx=0.5)


tune!(suite);
results = run(suite, verbose=true, seconds=100)


println("Writing to file...")
open("./paper/timing-results.txt", "w") do file

    println(file, "---")
    println(file, "Reflectance conversion for Dye_1-6")
    println(file, "---")
    println(file, mean(results["reflectance-conversion"]["Dye_1"]))
    println(file, "n-scanlines: ", read_envi_header(f1.bilhdr)["lines"])
    println(file, "\n")

    println(file, "---")
    println(file, "Reflectance conversion for NoDye_1")
    println(file, "---")
    println(file, mean(results["reflectance-conversion"]["NoDye_1"]))
    println(file, "n-scanlines: ", read_envi_header(f2.bilhdr)["lines"])
    println(file, "\n")

    println(file, "---")
    println(file, "Reflectance conversion for Scotty_1-1")
    println(file, "---")
    println(file, mean(results["reflectance-conversion"]["Scotty_1"]))
    println(file, "n-scanlines: ", read_envi_header(f3.bilhdr)["lines"])
    println(file, "\n")


    println(file, "The following timings used a cube of length ", read_envi_header(f3.bilhdr)["lines"])

    println(file, "---")
    println(file, "Resampling results for Δx=0.5")
    println(file, "---")
    println(file, mean(results["resampling"]["Δx=0.5"]))
    println(file, "\n")

    println(file, "---")
    println(file, "Resampling results for Δx=0.4")
    println(file, "---")
    println(file, mean(results["resampling"]["Δx=0.4"]))
    println(file, "\n")

    println(file, "---")
    println(file, "Resampling results for Δx=0.3")
    println(file, "---")
    println(file, mean(results["resampling"]["Δx=0.3"]))
    println(file, "\n")

    println(file, "---")
    println(file, "Resampling results for Δx=0.2")
    println(file, "---")
    println(file, mean(results["resampling"]["Δx=0.2"]))
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
sizes = [5, 10, 20, 30, 40, 50]
timings = [
    mean(results["resampling"]["Δx=0.05"]).time / 1e9,
    mean(results["resampling"]["Δx=0.1"]).time / 1e9,
    mean(results["resampling"]["Δx=0.2"]).time / 1e9,
    mean(results["resampling"]["Δx=0.3"]).time / 1e9,
    mean(results["resampling"]["Δx=0.4"]).time / 1e9,
    mean(results["resampling"]["Δx=0.5"]).time / 1e9,
]

scanlines = [
    parse(Int, read_envi_header(f2.bilhdr)["lines"]),
    parse(Int, read_envi_header(f1.bilhdr)["lines"]),
    parse(Int, read_envi_header(f3.bilhdr)["lines"]),
]

stimings = [
    mean(results["reflectance-conversion"]["NoDye_1"]).time / 1e9,
    mean(results["reflectance-conversion"]["Dye_1"]).time / 1e9,
    mean(results["reflectance-conversion"]["Scotty_1"]).time / 1e9,
]


my_theme = mints_theme
my_theme.Axis.xticklabelsize=20
my_theme.Axis.yticklabelsize=20
my_theme.Axis.xlabelsize=22
my_theme.Axis.ylabelsize=22
my_theme.Axis.titlesize=25
set_theme!(my_theme)

fig = Figure();
ax = CairoMakie.Axis(fig[1,1], xlabel="grid resolution (cm)", ylabel="execution time (seconds)", title="Hyperspectral Image Reinterpolation")
line  = lines!(ax, sizes, timings, linewidth=4)
scatter  = scatter!(ax, sizes, timings; markersize=15)
fig

save("paper/figures/regrid-timing.png", fig)
save("paper/figures/regrid-timing.pdf", fig)
save("paper/figures/regrid-timing.eps", fig)
save("paper/figures/regrid-timing.svg", fig)



fig = Figure();
ax = CairoMakie.Axis(fig[1,1], xlabel="number of scanlines", ylabel="execution time (seconds)", title="Loading & Reflectance Conversion")
line  = lines!(ax, scanlines, stimings, linewidth=4)
scatter  = scatter!(ax, scanlines, stimings; markersize=15)
fig

save("paper/figures/reflectance-timing.png", fig)
save("paper/figures/reflectance-timing.pdf", fig)
save("paper/figures/reflectance-timing.eps", fig)
save("paper/figures/reflectance-timing.svg", fig)


