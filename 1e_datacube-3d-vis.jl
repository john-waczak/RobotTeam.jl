using RobotTeam
using HDF5
using BenchmarkTools
using GLMakie, GeometryBasics
using MintsMakieRecipes
using ProgressMeter

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





# 1. create visualization of un-georectified HSI
size_in_inches = (3, 3)
dpi = 400
size_in_pixels = size_in_inches .* dpi

fig = vis_cube(hsi; cmap=:jet, offset=0.1, ibounds=(250,1600), resolution=size_in_pixels)

save("./paper/figures/demo-cube.png", fig)


# 2. create visualization of georectified HSI

# now let's visualize the georectified cube
xs, ys, isnorth, zone, Longitudes, Latitudes, IsInbounds, varnames, printnames, λs, Data_μ = resample_datacube_fast(hsi; Δx=Δx)


fig = vis_rectified_cube(Data_μ, xs, ys, IsInbounds, λs, Δx; offset=0.0, ibounds=(100, length(xs)), resolution=size_in_pixels, colorbar=true)

colgap!(fig.layout, -150)

fig

save("./paper/demo-rectified-cube.png", fig)



