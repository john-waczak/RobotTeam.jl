using RobotTeam
using HDF5
using BenchmarkTools
using CairoMakie
using MintsMakieRecipes

set_theme!(mints_theme)

my_theme = mints_theme
my_theme.Axis.xticklabelsize=17
my_theme.Axis.yticklabelsize=17
my_theme.Axis.xlabelsize=20
my_theme.Axis.ylabelsize=20
my_theme.Axis.titlesize=22
set_theme!(my_theme)


include("utils/vis_tools.jl")
include("utils/config.jl")




basepath


# create a plot of the HSI path on a satellite background tile
w= -97.717472
n= 33.703572
s= 33.700797
e= -97.712413


satmap = get_background_satmap(w,e,s,n)


h5path = joinpath(outpath, "12-09", "Dye_1", "Dye_1-6.h5")
@assert isfile(h5path)


h5 = h5open(h5path, "r")

Latitudes = read(h5["data-Δx_0.1/Latitudes"])
Longitudes = read(h5["data-Δx_0.1/Longitudes"])
Data_μ = read(h5["data-Δx_0.1/Data_μ"])
varnames = read(h5["data-Δx_0.1/varnames"])
printnames = read(h5["data-Δx_0.1/printnames"])
λs = read(h5["data-Δx_0.1/λs"])


indices_metrics = findfirst(x-> x=="mNDWI", varnames):length(varnames)

spec_R = Data_μ[1:]

size_in_inches = (4, 3)
dpi = 300
size_in_pixels = size_in_inches .* dpi

fig = Figure(;resolution=size_in_pixels)
ax = CairoMakie.Axis(fig[1,1], xlabel="Longitude", ylabel="Latitude")
hm = heatmap!(ax, satmap.w..satmap.e, satmap.s..satmap.n, satmap.img)

rgb = getRGB(h5)
hm_rgb = heatmap!(ax, Longitudes, Latitudes, rgb)

xlims!(ax, -97.717, -97.715)
ylims!(ax, 33.7015, 33.7025)

save("./paper/figures/georectified-hsi-plume.png", fig)

fig

close(h5)

# plot radiance, irradiance, and reflectance for single HSI (12-09 Dye_1-6) for a slice from the plume
fpath = "/media/teamlary/LabData/RobotTeam/raw/hsi/12-09/Dye_1-6"
readdir(fpath)
@assert ispath(fpath)

rad_dat = zeros(UInt16, length(λs))
ref_dat = zeros(Float64, length(λs))
let
    dat, h, p = read_envi_file(joinpath(fpath, "Dye_1_Pika_XC2_6-radiance.bil"), joinpath(fpath, "Dye_1_Pika_XC2_6-radiance.bil.hdr"))
    rad_dat .= dat[:, 100, 100]
end

fnames = joinpath.(fpath, readdir(fpath))

let
    hsi = HyperspectralImage(
        fnames[1],
        fnames[2],
        fnames[4],
        fnames[3],
        fnames[5],
        fnames[6];
        isflipped=true
    );

    ref_dat .= hsi.Reflectance[1:length(λs), 100, 100]
end


irrad_dat = zeros(Float64, 2048)
λs_irrad = zeros(Float64, 2048)
let
    dat, h, p = read_envi_file(fnames[5], fnames[6])
    irrad_dat .= vec(dat)
    λs_irrad .= p["wavelengths"]
end


fig = Figure()
ax = CairoMakie.Axis(fig[1,1], xlabel="λ (nm)", ylabel="Radiance", yticksvisible=false, ygridvisible=false, yminorgridvisible=false, yticklabelsvisible=false)
ax2 = CairoMakie.Axis(fig[1,1], yaxisposition=:right, ylabel="Downwelling Irradiance", yticksvisible=false, ygridvisible=false, yminorgridvisible=false, yticklabelsvisible=false)
linkxaxes!(ax, ax2)

ylims!(ax, 0, nothing)
xlims!(ax, min(λs[1], λs_irrad[1]), max(λs[end], λs_irrad[end]))
ylims!(ax2, 0, nothing)
xlims!(ax2, min(λs[1], λs_irrad[1]), max(λs[end], λs_irrad[end]))


lin = lines!(ax, λs, rad_dat, color=:gray, linewidth=3, label="Radiance")

lin2 = lines!(ax2, λs_irrad, irrad_dat, color=:orange, linewidth=3, label="Downwelling Irradiance")

leg = Legend(fig[1,1], [lin, lin2], ["Radiance", "Downwelling Irradiance"], orientation=:horizontal,  tellheight = false, tellwidth = false, halign=:right, valign=:top, margin = (10, 30, 10, 10),)

fig

save("./paper/figures/radiance-sample.png", fig)
save("./paper/figures/radiance-sample.svg", fig)
save("./paper/figures/radiance-sample.eps", fig)
save("./paper/figures/radiance-sample.pdf", fig)


λ_blue = 495.0
λ_green = 530.0
λ_red = 650.0

cg = cgrad([:purple, :blue, :green, :red, :maroon], [0, λ_blue/1000, λ_green/1000, λ_red/1000  ,1])



fig = Figure()
ax = CairoMakie.Axis(fig[1,1], xlabel="λ (nm)", ylabel="Reflectance", yticksvisible=false, ygridvisible=false, yminorgridvisible=false, yticklabelsvisible=false)

colors = [cg[λ/maximum(λs)] for λ ∈ λs]
lin = lines!(ax, λs, ref_dat, linewidth=3, color=colors)

ylims!(ax, 0, nothing)
xlims!(ax, λs[1], λs[end])

save("./paper/figures/reflectance-sample.png", fig)
save("./paper/figures/reflectance-sample.svg", fig)
save("./paper/figures/reflectance-sample.eps", fig)
save("./paper/figures/reflectance-sample.pdf", fig)
fig

# plot single HSI on map (1209 Dye_1-6)

# plot spectral indices filtered to water-only pixels


