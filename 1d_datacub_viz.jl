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


indices_metrics = findfirst(x-> x=="mNDWI", varnames):length(varnames)


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

# plot radiance, irradiance, and reflectance for single HSI (12-09 Dye_1-6) for a slice from the plume


# plot single HSI on map (1209 Dye_1-6)

# plot spectral indices filtered to water-only pixels


