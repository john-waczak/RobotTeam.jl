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


# fetch georectified image
hsi = HyperspectralImage(paths...; isflipped=true)  # 4.104 s


# test color image:
img = getRGB(hsi)
cmk.image(img)


viz(mesh, color=vcat(img...))


# create a plot of the HSI path on a satellite background tile
w= -97.717472
n= 33.703572
s= 33.700797
e= -97.712413


satmap = get_background_satmap(w,e,s,n)
mesh, ∂mesh = build_mesh(hsi.Longitudes .- satmap.w, hsi.Latitudes .- satmap.s)
img = getRGB(hsi)

f = cmk.Figure(resolution=(800,600))
ax = cmk.Axis(f[1,1], xlabel="longitude", ylabel="latitude")
cmk.heatmap!(ax, 0.0..(satmap.e-satmap.w), 0.0..(satmap.n-satmap.s), satmap.img)
cmk.xlims!(ax, -97.7168-satmap.w, -97.7125-satmap.w)
cmk.ylims!(ax, 33.70075-satmap.s, 33.7035-satmap.s)

#viz!(ax, mesh, color=vcat(img...))
viz!(ax, mesh, color=vcat(hsi.SolarAzimuth...))

f

# now we want to resample
