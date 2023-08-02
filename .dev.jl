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

basepath = "/media/jwaczak/Data/robotteam-data/sample-hsi/Scotty_1-2"

outpath = "/media/jwaczak/Data/robotteam-data/h5"
if !ispath(outpath)
    mkpath(outpath)
end

paths = [
    joinpath(basepath, "Scotty_1_Pika_XC2_2-radiance.bil"),
    joinpath(basepath, "Scotty_1_Pika_XC2_2-radiance.bil.hdr"),
    joinpath(basepath, "Scotty_1_Pika_XC2_2.lcf"),
    joinpath(basepath, "Scotty_1_Pika_XC2_2.bil.times" ),
    joinpath(basepath, "Scotty_1_downwelling_2_pre.spec"),
    joinpath(basepath, "Scotty_1_downwelling_2_pre.spec.hdr")
]
@assert all(ispath.(paths))

h5path = joinpath(outpath, "test_hsi.h5")

# return opened h5 file
h5 = envi_to_hdf5(paths..., h5path)

generateReflectance!(h5)
generateDerivedMetrics!(h5)


#generateCoords!(h5, θ_view=30.8, z_ground=292.0, isflipped=true)
generateCoords!(h5, θ_view=30.8, z_ground=292.0, isflipped=true)
mesh, ∂mesh = build_mesh(read(h5["raw/georectified/longitude"]), read(h5["raw/georectified/latitude"]));

λs = read(h5["raw/reflectance/wavelengths"])
idx_r = argmin(abs.(λs .- 630.0))
idx_g = argmin(abs.(λs .- 532.0))
idx_b = argmin(abs.(λs .- 465.0))

img = permutedims(read(h5["raw/reflectance/reflectance"])[[idx_r, idx_g, idx_b],:,:], (1,3,2))

imgp = process_image(img)


cmk.image(rotr90(colorview(RGB, imgp)),
          axis = (aspect = cmk.DataAspect(), yreversed=false,),
          )



# create a plot of the HSI path on a satellite background tile
lon = read(h5["raw/lcf/longitude"])
lat = read(h5["raw/lcf/latitude"])
lons = read(h5["raw/georectified/longitude"])
lats = read(h5["raw/georectified/latitude"])

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

hsi = viz!(ax, mesh, color=vcat(colorview(RGB, imgp)...))
∂hsi = viz!(ax, ∂mesh, color=:red)
flightpath = cmk.lines!(ax, lon, lat, color=:white)

f

save("test_georectification.png", f)



close(h5)
