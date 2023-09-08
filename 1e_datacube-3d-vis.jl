using RobotTeam
using HDF5
using BenchmarkTools
using GLMakie, GeometryBasics
using MintsMakieRecipes
using ProgressMeter
using Statistics

set_theme!(mints_theme)
my_theme = mints_theme
my_theme.Axis.xticklabelsize = 20
my_theme.Axis.yticklabelsize = 20
my_theme.Axis.xlabelsize = 22
my_theme.Axis.ylabelsize = 22
my_theme.Axis.titlesize = 25
my_theme.Colorbar.ticklabelsize=20
my_theme.Colorbar.labelsize=22
set_theme!(my_theme)



include("utils/vis_tools.jl")
include("utils/config.jl")



basepath = "/home/teamlary/gitrepos/utd/RobotTeam.jl/data/raw"
f = get_raw_file_list(get_bil_files(basepath, "Dye_1")[1])


# read in the HSI
hsi = HyperspectralImage(f.bilpath, f.bilhdr, f.lcfpath, f.timespath, f.specpath, f.spechdr; isflipped=true);

img = getRGB(hsi)
fig, ax, hm = heatmap(img)
size(img)


# 1. create visualization of un-georectified HSI
size_in_inches = (3, 3)
dpi = 600
size_in_pixels = size_in_inches .* dpi

fig = vis_cube(hsi; cmap=:jet, offset=0.1, ibounds=(250,1600), resolution=size_in_pixels)

save("./paper/figures/demo-cube.png", fig)


# 2. create visualization of georectified HSI

# now let's visualize the georectified cube
xs, ys, isnorth, zone, Longitudes, Latitudes, IsInbounds, varnames, printnames, λs, Data_μ = resample_datacube_fast(hsi; Δx=Δx)


size_in_inches = (3, 3)
dpi = 300
size_in_pixels = size_in_inches .* dpi



fig = vis_rectified_cube(Data_μ, xs, ys, IsInbounds, λs, Δx; offset=0.0, ibounds=(100, length(xs)), resolution=size_in_pixels, colorbar=true)

colgap!(fig.layout, -100)

fig

save("./paper/figures/demo-rectified-cube.png", fig)



# plot reflectance of slice through the plume
using CairoMakie

set_theme!(figure_padding = 30)

data = log10.(Data_μ[1:462,:,:])
idx_not_nan_or_inf = findall(.!(isnan.(data)) .&& .!(isinf.(data)))
Rmin = quantile(data[idx_not_nan_or_inf], 0.1)
Rmax = quantile(data[idx_not_nan_or_inf], 0.99)



ij_pixels = findall(IsInbounds)
rgb = getRGB(Data_μ, λs, ij_pixels)

fig1, ax, hm = heatmap(rgb)
scatter!(ax, [80], [225], marker=:circle, color=:cyan, markersize=10)
scatter!(ax, [400], [320], marker=:circle, color=:cyan, markersize=10)
scatter!(ax, [185], [45], marker=:circle, color=:cyan, markersize=10)
scatter!(ax, [400], [80], marker=:circle, color=:cyan, markersize=10)
fig1

hidedecorations!(ax)
hidespines!(ax)

fig1

save("./paper/figures/hsi-points-of-interest.png", fig1)
save("./paper/figures/hsi-points-of-interest.eps", fig1)
save("./paper/figures/hsi-points-of-interest.pdf", fig1)


idx_plume = CartesianIndex(80, 225)
idx_water = CartesianIndex(400, 320)
idx_algae = CartesianIndex(186, 45)
idx_grass = CartesianIndex(400, 80)


fig = Figure();
ax = CairoMakie.Axis(fig[1,1], xlabel="λ (nm)", ylabel="Reflectance");
ls_plume = lines!(ax, λs, Data_μ[1:462, idx_plume], linewidth=3, color=:red)
ls_water = lines!(ax, λs, Data_μ[1:462, idx_water], linewidth=3, color=:blue)
ls_algae = lines!(ax, λs, Data_μ[1:462, idx_algae], linewidth=3, color=:green)
ls_grass = lines!(ax, λs, Data_μ[1:462, idx_grass], linewidth=3, color=:brown)
leg = axislegend(
    ax,
    [ls_plume,ls_water, ls_algae, ls_grass],
    ["Rhodamine Plume", "Open Water", "Algae", "Grass"],
    position=:lt,
)
xlims!(λs[1], λs[end])
fig

save("./paper/figures/reflectance-samples.png", fig)
save("./paper/figures/reflectance-samples.eps", fig)
save("./paper/figures/reflectance-samples.svg", fig)
save("./paper/figures/reflectance-samples.pdf", fig)


# plume
fig2 = Figure();
ax = CairoMakie.Axis(fig2[1, 1], xlabel="λ (nm)", ylabel="Reflectance");
ls = lines!(ax, λs, Data_μ[1:462, idx_plume], color=data[1:462, idx_plume], colormap=:jet, colorrange=(Rmin, Rmax), linewidth=3.5)
xlims!(λs[1], λs[end])
fig2

save("./paper/figures/reflectance-plume.png", fig2)
save("./paper/figures/reflectance-plume.eps", fig2)
save("./paper/figures/reflectance-plume.svg", fig2)
save("./paper/figures/reflectance-plume.pdf", fig2)


# water
fig2 = Figure();
ax = CairoMakie.Axis(fig2[1, 1], xlabel="λ (nm)", ylabel="Reflectance");
ls = lines!(ax, λs, Data_μ[1:462, idx_water], color=data[1:462, idx_water], colormap=:jet, colorrange=(Rmin, Rmax), linewidth=3.5)
xlims!(λs[1], λs[end])
fig2

save("./paper/figures/reflectance-water.png", fig2)
save("./paper/figures/reflectance-water.eps", fig2)
save("./paper/figures/reflectance-water.svg", fig2)
save("./paper/figures/reflectance-water.pdf", fig2)


# algae
# see: https://www.researchgate.net/figure/Comparison-of-reflectance-spectra-for-algal-taxa-representing-four-different-divisions_fig11_226194233

fig2 = Figure();
ax = CairoMakie.Axis(fig2[1, 1], xlabel="λ (nm)", ylabel="Reflectance");
ls = lines!(ax, λs, Data_μ[1:462, idx_algae], color=data[1:462, idx_algae], colormap=:jet, colorrange=(Rmin, Rmax), linewidth=3.5)
xlims!(λs[1], λs[end])
fig2
# xlims!(400, 800)
# ylims!(0, 0.03)

# fig2

save("./paper/figures/reflectance-algae.png", fig2)
save("./paper/figures/reflectance-algae.eps", fig2)
save("./paper/figures/reflectance-algae.svg", fig2)
save("./paper/figures/reflectance-algae.pdf", fig2)



# grass
fig2 = Figure();
ax = CairoMakie.Axis(fig2[1, 1], xlabel="λ (nm)", ylabel="Reflectance");
ls = lines!(ax, λs, Data_μ[1:462, idx_grass], color=data[1:462, idx_grass], colormap=:jet, colorrange=(Rmin, Rmax), linewidth=3.5)
xlims!(λs[1], λs[end])
fig2

save("./paper/figures/reflectance-grass.png", fig2)
save("./paper/figures/reflectance-grass.eps", fig2)
save("./paper/figures/reflectance-grass.svg", fig2)
save("./paper/figures/reflectance-grass.pdf", fig2)

















