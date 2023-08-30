using RobotTeam
using HDF5
using BenchmarkTools
using GLMakie, GeometryBasics
using MintsMakieRecipes
using ProgressMeter
using Statistics

set_theme!(mints_theme)

include("utils/vis_tools.jl")
include("utils/config.jl")



basepath = "/home/teamlary/gitrepos/utd/RobotTeam.jl/data/raw"
f = get_raw_file_list(get_bil_files(basepath, "Dye_1")[1])


# read in the HSI
hsi = HyperspectralImage(f.bilpath, f.bilhdr, f.lcfpath, f.timespath, f.specpath, f.spechdr; isflipped=true);

img = getRGB(hsi)
fig, ax, hm = heatmap(img)
size(img)



my_theme = mints_theme
my_theme.Axis.xticklabelsize=17
my_theme.Axis.yticklabelsize=17
my_theme.Axis.xlabelsize=20
my_theme.Axis.ylabelsize=20
my_theme.Axis.titlesize=22
set_theme!(my_theme)



# 1. create visualization of un-georectified HSI
size_in_inches = (3, 3)
dpi = 600
size_in_pixels = size_in_inches .* dpi

fig = vis_cube(hsi; cmap=:jet, offset=0.1, ibounds=(250,1600), resolution=size_in_pixels)

save("./paper/figures/demo-cube.png", fig)


# 2. create visualization of georectified HSI

# now let's visualize the georectified cube
xs, ys, isnorth, zone, Longitudes, Latitudes, IsInbounds, varnames, printnames, λs, Data_μ = resample_datacube_fast(hsi; Δx=Δx)


fig = vis_rectified_cube(Data_μ, xs, ys, IsInbounds, λs, Δx; offset=0.0, ibounds=(100, length(xs)), resolution=size_in_pixels, colorbar=true)

colgap!(fig.layout, -300)

fig

save("./paper/figures/demo-rectified-cube.png", fig)


ax = current_axis(fig)

arrows!([(length(xs)-250)*Δx],[(length(ys))*Δx],[Δx*(-60)],[0],[0],[50*Δx]; arrowsize=1)


save("./paper/figures/rectified-cube-w-arrow.png", fig)

fig

# plot reflectance of slice through the plume
using CairoMakie

data = log10.(Data_μ[1:462,:,:])
idx_not_nan_or_inf = findall(.!(isnan.(data)) .&& .!(isinf.(data)))
Rmin = quantile(data[idx_not_nan_or_inf], 0.1)
Rmax = quantile(data[idx_not_nan_or_inf], 0.99)

fig2 = Figure();
ax = GLMakie.Axis(fig2[1,1], xlabel="λ (nm)", ylabel="Reflectance");
ls = lines!(ax, λs, Data_μ[1:462,80, end-150], color= data[1:462,80, end-150], colormap=:jet, colorrange=(Rmin, Rmax), linewidth=3.5)
xlims!(λs[1], λs[end])
fig2

save("./paper/figures/reflectance-at-arrow.png", fig2)
save("./paper/figures/reflectance-at-arrow.eps", fig2)
save("./paper/figures/reflectance-at-arrow.svg", fig2)
save("./paper/figures/reflectance-at-arrow.pdf", fig2)


# ij_pixels = findall(IsInbounds)
# Ref_img = getRGB(Data_μ, λs, ij_pixels)
# fig2, ax, hm = heatmap(xs[1]..xs[end], ys[1]..ys[end], Ref_img)
# scatter!(ax, xs[80],ys[end-150],marker=:circle, color=:cyan, markersize=30)
# fig2



