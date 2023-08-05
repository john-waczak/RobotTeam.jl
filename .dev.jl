using RobotTeam
using HDF5
using BenchmarkTools
using DelimitedFiles
using Dates

import CairoMakie as cmk
using IntervalSets
using MintsMakieRecipes
cmk.set_theme!(mints_theme)

using LoopVectorization
using LinearAlgebra
using StaticArrays

using Meshes, MeshViz

using Images, ImageContrastAdjustment

include("utils/vis_tools.jl")



# 1. set up test paths

basepath = "/media/jwaczak/Data/robotteam-data/sample-hsi/Scotty_1-1"

outpath = "/media/jwaczak/Data/robotteam-data/h5"
if !ispath(outpath)
    mkpath(outpath)
end

paths = [
    joinpath(basepath, "Scotty_1_Pika_XC2_1-radiance.bil"),
    joinpath(basepath, "Scotty_1_Pika_XC2_1-radiance.bil.hdr"),
    joinpath(basepath, "Scotty_1_Pika_XC2_1.lcf"),
    joinpath(basepath, "Scotty_1_Pika_XC2_1.bil.times" ),
    joinpath(basepath, "Scotty_1_downwelling_1_pre.spec"),
    joinpath(basepath, "Scotty_1_downwelling_1_pre.spec.hdr")
]
@assert all(ispath.(paths))

h5path = joinpath(outpath, "test_hsi.h5")


# set up data readers for individual HSI


h5 = h5open(h5path)
size(h5["raw/radiance/radiance"])
close(h5)


typeof(zeros(2,2,2))



bilpath = paths[1]
bilhdr = paths[2]
lcfpath = paths[3]
timespath = paths[4]
specpath = paths[5]
spechdr = paths[6]


# 0. read header
h = read_envi_header(bilhdr)
p = get_envi_params(h)


# 1. preallocate HSI struct
hsi = HyperspectralImage(p["nbands"], p["ncols"], p["nrows"])
# @benchmark HyperspectralImage(size(img)...)

# 2. load radiance
img, h, p = read_envi_file(bilpath, bilhdr)
# @benchmark read_envi_file(bilpath, bilhdr)  # 444.6 ms
year = parse(Int, split(split(p["timestamp"], "-")[1], "/")[end])

hsi.Radiance .=  img;
hsi.λs .= p["wavelengths"]
 
# 2. compute Reflectance
generateReflectance!(hsi, specpath, spechdr)

# 3. generate flight data
fdata = FlightData(lcfpath, timespath)

# 4. georectify coordinates
generateCoords!(hsi, fdata;isflipped=true)


extrema(hsi.X[1,:,:])
extrema(hsi.X[2,:,:])


mesh, ∂mesh = build_mesh(hsi.X[1,:,:], hsi.X[2,:,:]);

f = cmk.Figure(resolution=(1600,800))
ax = cmk.Axis(f[1,1])
viz!(ax, ∂mesh)

ax2 = cmk.Axis(f[1,2])
mesh2, ∂mesh2 = build_mesh(hsi.X[1,:,:] .- mean(hsi.X[1,:,:]), hsi.X[2,:,:] .- mean(hsi.X[2,:,:]));
viz!(ax2, ∂mesh2)

f



# test color image:
img = getRGB(hsi)
cmk.image(img)


viz(mesh2, color=vcat(img...))





# mesh, ∂mesh = build_mesh(hsi.X[1,:,:], hsi.X[2,:,:]);

# viz(∂mesh)


# create a plot of the HSI path on a satellite background tile
w= -97.717472
n= 33.703572
s= 33.700797
e= -97.712413


satmap = get_background_satmap(w,e,s,n)


f = cmk.Figure(resolution=(800,600))

ax = cmk.Axis(f[1,1], xlabel="longitude", ylabel="latitude")

cmk.heatmap!(ax, satmap.w..satmap.e, satmap.s..satmap.n, satmap.img)
cmk.xlims!(ax, -97.7168, -97.7125)
cmk.ylims!(ax, 33.70075, 33.7035)

f

hsi_plot = viz!(ax, mesh, color=vcat(img...))
∂hsi_plot = viz!(ax, ∂mesh, color=:red)
# save("test_georectification.png", f)

f


