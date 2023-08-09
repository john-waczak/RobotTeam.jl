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


hsi = HyperspectralImage(paths...; isflipped=true)  # 4.104 s


# create vector of labels
varnames = []
printnames = []
for i ∈ 1:length(hsi.λs)
    idx = lpad(i, 3, "0")
    push!(varnames, "R_$(idx)")
    push!(printnames, "Reflectance Band $(idx)")
end

push!(varnames, "roll")
push!(printnames, "Roll")

push!(varnames, "pitch")
push!(printnames, "Pitch")

push!(varnames, "heading")
push!(printnames, "Heading")

push!(varnames, "camera_angle")
push!(printnames, "Camera Angle")

push!(varnames, "solar_azimuth")
push!(printnames, "Solar Azimuth")

push!(varnames, "solar_elevation")
push!(printnames, "Solar Elevation")

push!(varnames, "solar_zenith")
push!(printnames, "Solar Zenith")

push!(varnames, "mNDWI")
push!(printnames, "mNDWI")

push!(varnames, "NDVI")
push!(printnames, "NDVI")

push!(varnames, "SR")
push!(printnames, "SR")

push!(varnames, "EVI")
push!(printnames, "EVI")

push!(varnames, "AVRI")
push!(printnames, "AVRI")

push!(varnames, "NDVI_705")
push!(printnames, "NDVI_705")

push!(varnames, "MSR_705")
push!(printnames, "MSR_705")

push!(varnames, "MNDVI")
push!(printnames, "MNDVI")

push!(varnames, "VOG1")
push!(printnames, "VOG1")

push!(varnames, "VOG2")
push!(printnames, "VOG2")

push!(varnames, "VOG3")
push!(printnames, "VOG3")

push!(varnames, "PRI")
push!(printnames, "PRI")

push!(varnames, "SIPI")
push!(printnames, "SIPI")

push!(varnames, "PSRI")
push!(printnames, "PSRI")

push!(varnames, "CRI1")
push!(printnames, "CRI1")

push!(varnames, "CRI2")
push!(printnames, "CRI2")

push!(varnames, "ARI1")
push!(printnames, "ARI1")

push!(varnames, "ARI2")
push!(printnames, "ARI2")

push!(varnames, "WBI")
push!(printnames, "WBI")

push!(varnames, "MCRI")
push!(printnames, "MCRI")

push!(varnames, "TCARI")
push!(printnames, "TCARI")





# mesh, ∂mesh = build_mesh(hsi.X[1,:,:], hsi.X[2,:,:]);

# f = cmk.Figure(resolution=(1600,800))
# ax = cmk.Axis(f[1,1])
# viz!(ax, ∂mesh)

# ax2 = cmk.Axis(f[1,2])
# mesh2, ∂mesh2 = build_mesh(hsi.X[1,:,:] .- mean(hsi.X[1,:,:]), hsi.X[2,:,:] .- mean(hsi.X[2,:,:]));
# viz!(ax2, ∂mesh2)

# f


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

f = cmk.Figure(resolution=(800,600))

ax = cmk.Axis(f[1,1], xlabel="longitude", ylabel="latitude")

cmk.heatmap!(ax, 0.0..(satmap.e-satmap.w), 0.0..(satmap.n-satmap.s), satmap.img)
f
cmk.xlims!(ax, -97.7168-satmap.w, -97.7125-satmap.w)
cmk.ylims!(ax, 33.70075-satmap.s, 33.7035-satmap.s)

mesh, ∂mesh = build_mesh(hsi.Longitudes .- satmap.w, hsi.Latitudes .- satmap.s)

img = getRGB(hsi)
viz!(mesh, color=vcat(img...))

f




