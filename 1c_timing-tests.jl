using RobotTeam
using HDF5
using BenchmarkTools
using CairoMakie
using MintsMakieRecipes

set_theme!(mints_theme)

update_theme!(
    figure_padding=30,
    Axis=(
        xticklabelsize=20,
        yticklabelsize=20,
        xlabelsize=22,
        ylabelsize=22,
        titlesize=25,
    ),
    Colorbar=(
        ticklabelsize=20,
        labelsize=22
    )
)


include("utils/vis_tools.jl")
include("utils/config.jl")


basepath = "/home/teamlary/gitrepos/utd/RobotTeam.jl/data/raw"
f1 = get_raw_file_list(get_bil_files(basepath, "NoDye_1")[1])
f2 = get_raw_file_list(get_bil_files(basepath, "NoDye_1")[2])
f3 = get_raw_file_list(get_bil_files(basepath, "NoDye_1")[3])

f = f3


sizes = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
scanlines = [
    parse(Int, read_envi_header(f1.bilhdr)["lines"]),
    parse(Int, read_envi_header(f2.bilhdr)["lines"]),
    parse(Int, read_envi_header(f3.bilhdr)["lines"]),
]


georeferencing_timings = []
resampling_timings = []
reflectance_timings = []

# benchmark the loading and georeferencing timings
let
    t = @benchmark HyperspectralImage($f1.bilpath, $f1.bilhdr, $f1.lcfpath, $f1.timespath; isflipped=true)
    push!(georeferencing_timings, mean(t.times)/1e9)
end

let
    t = @benchmark HyperspectralImage($f2.bilpath, $f2.bilhdr, $f2.lcfpath, $f2.timespath; isflipped=true)
    push!(georeferencing_timings, mean(t.times)/1e9)
end

let
    t = @benchmark HyperspectralImage($f3.bilpath, $f3.bilhdr, $f3.lcfpath, $f3.timespath; isflipped=true)
    push!(georeferencing_timings, mean(t.times)/1e9)
end



# benchmark resampling times
let
    hsi = HyperspectralImage(f.bilpath, f.bilhdr, f.lcfpath, f.timespath; isflipped=true);

    for i in 1:length(sizes)
        size = sizes[i]/100
        t = @benchmark resample_datacube_fast($hsi; Δx=$size)
        push!(resampling_timings, mean(t.times)/1e9)
    end
end



# benchmark reflectance timings
let
    hsi1 = HyperspectralImage(f1.bilpath, f1.bilhdr, f1.lcfpath, f1.timespath; isflipped=true);
    _, _, _, _, _, _, _, varnames, printnames, λs, Data1 = resample_datacube_fast(hsi1);

    size(Data1)
    # (469, 595, 396)
    # 936 ms
    t = @benchmark generateReflectance!($Data1, $f1.specpath, $f1.spechdr, $λs)
    push!(reflectance_timings, mean(t.times)/1e9)
end

let
    hsi2 = HyperspectralImage(f2.bilpath, f2.bilhdr, f2.lcfpath, f2.timespath; isflipped=true);
    _, _, _, _, _, _, _, varnames, printnames, λs, Data2 = resample_datacube_fast(hsi2)

    size(Data2)
    # 114 ms
    # (469, 620, 389)
    t = @benchmark generateReflectance!($Data2, $f2.specpath, $f2.spechdr, $λs)
    push!(reflectance_timings, mean(t.times)/1e9)
end

let
    hsi3 = HyperspectralImage(f3.bilpath, f3.bilhdr, f3.lcfpath, f3.timespath; isflipped=true);
    _, _, _, _, _, _, _, varnames, printnames, λs, Data3 = resample_datacube_fast(hsi3)

    size(Data3)
    # (469, 487, 947)
    # 161 ms
    t = @benchmark generateReflectance!($Data3, $f3.specpath, $f3.spechdr, $λs)
    push!(reflectance_timings, mean(t.times)/1e9)
end


resampling_timings = Float64.(resampling_timings)
georeferencing_timings = Float64.(georeferencing_timings)
reflectance_timings = Float64.(reflectance_timings)

fig = Figure();
ax = CairoMakie.Axis(fig[1,1], xlabel="Grid Resolution (cm)", ylabel="Execution Time (seconds)", title="Hyperspectral Image Resampling")
line  = lines!(ax, sizes, resampling_timings, linewidth=4)
scatter  = scatter!(ax, sizes, resampling_timings; markersize=15)
xlims!(5, 50)
fig

save("paper/figures/regrid-timing.png", fig)
save("paper/figures/regrid-timing.pdf", fig)
save("paper/figures/regrid-timing.eps", fig)
save("paper/figures/regrid-timing.svg", fig)


size_in_inches = (5, 3)
dpi = 300
size_in_pixels = size_in_inches .* dpi

scanlines
reflectance_timings

fig = Figure(resolution=size_in_pixels);
ax1 = CairoMakie.Axis(fig[1,1], xlabel="number of scanlines", ylabel="Execution Time (seconds)", title="Loading and Georeferencing")
ax2 = CairoMakie.Axis(fig[1,2], xlabel="number of scanlines", ylabel="Execution Time (seconds)", title="Radiance to Reflectance Conversion", yaxisposition=:right)

line1  = lines!(ax1, scanlines, georeferencing_timings; linewidth=4)
scatter1  = scatter!(ax1, scanlines, georeferencing_timings; markersize=15)

line2  = lines!(ax2, scanlines, reflectance_timings; linewidth=4, color=mints_colors[2])
scatter2 = scatter!(ax2, scanlines, reflectance_timings; markersize=15, color=mints_colors[2])

fig


save("paper/figures/reflectance-timing.png", fig)
save("paper/figures/reflectance-timing.pdf", fig)
save("paper/figures/reflectance-timing.eps", fig)
save("paper/figures/reflectance-timing.svg", fig)


println("Writing to file...")
open("./paper/timing-results.txt", "w") do file


    println(file, "\n--------------------------\n")
    println(file, "Resampling timings for a cube of length ", read_envi_header(f3.bilhdr)["lines"], ":")
    println(file, "\n--------------------------\n")

    for i ∈ 1:length(sizes)
        println(file, "---")
        println(file, "Resampling results for Δx=$(sizes[i])")
        println(file, "---")
        println(file, resampling_timings[i])

    end

    println(file, "\n--------------------------\n")
    println(file, "Georeferencing and Loading:")
    println(file, "\n--------------------------\n")

    for i ∈ 1:length(scanlines)
        println(file, "---")
        println(file, "Loading and Georeferencing time for $(scanlines[i]) scanlines")
        println(file, "---")
        println(file, georeferencing_timings[i])
        println(file, "\n")
    end

    println(file, "\n--------------------------\n")
    println(file, "Reflectance conversion:")
    println(file, "\n--------------------------\n")

    for i ∈ 1:length(scanlines)
        println(file, "---")
        println(file, "Reflectance conversion time for $(scanlines[i]) scanlines")
        println(file, "---")
        println(file, reflectance_timings[i])
        println(file, "\n")
    end
end


