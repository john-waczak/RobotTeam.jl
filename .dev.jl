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

using Meshes




# 1. set up test paths

basepath = "/media/jwaczak/Data/robotteam-data/sample-hsi/NoDye_1-1"
outpath = "/media/jwaczak/Data/robotteam-data/h5"
if !ispath(outpath)
    mkpath(outpath)
end


@assert ispath(basepath)

bil_path = joinpath(basepath, "NoDye_1_Pika_XC2_1-radiance.bil")
@assert ispath(bil_path)

bil_hdr = joinpath(basepath, "NoDye_1_Pika_XC2_1-radiance.bil.hdr")
@assert ispath(bil_hdr)

bil_times = joinpath(basepath, "NoDye_1_Pika_XC2_1.bil.times" )
@assert ispath(bil_times)

bil_lcf = joinpath(basepath, "NoDye_1_Pika_XC2_1.lcf")
@assert ispath(bil_lcf)

spec_path = joinpath(basepath, "NoDye_1_downwelling_1_pre.spec")
@assert ispath(spec_path)

spec_hdr = joinpath(basepath, "NoDye_1_downwelling_1_pre.spec.hdr")
@assert ispath(spec_hdr)


h5path = joinpath(outpath, "NoDye_1-1.h5")


# return opened h5 file
h5 = envi_to_hdf5(bil_path, bil_hdr, bil_lcf, bil_times, spec_path, spec_hdr, h5path)


generateReflectance!(h5)
generateDerivedMetrics!(h5)
generateCoords!(h5, θ_view=30.8, z_ground=292.0, isflipped=true)


# update generateCoords to convert utmz back into lat/lon

h5
close(h5)



# create a plot of the HSI path on a satellite background tile
lon = read(h5["raw/lcf/longitude"])
lat = read(h5["raw/lcf/latitude"])
times = read(h5["raw/lcf/times"])

roll = read(h5["raw/lcf/roll"])
pitch = read(h5["raw/lcf/pitch"])
yaw = read(h5["raw/lcf/heading"])

w= -97.717472
n= 33.703572
s= 33.700797
e= -97.712413


basemap = get_basemap(w,e,s,n)


f = cmk.Figure(resolution=(1200,500))

gleft = f[1,1] = cmk.GridLayout()
ax = cmk.Axis(gleft[1,1], xlabel="longitude", ylabel="latitude")
cmk.heatmap!(ax, basemap.w..basemap.e, basemap.s..basemap.n, basemap.map)
cmk.xlims!(ax, -97.7165, -97.7125)
cmk.ylims!(ax, 33.7011, 33.7032)

ls = cmk.lines!(ax, lon, lat, color=times)

lon
lat

cb = cmk.Colorbar(gleft[1,2], ls, label="time (s)")

f

gright = f[1,2] = cmk.GridLayout()

ax1 = cmk.Axis(gright[1,1], ylabel="roll (°)")
ax2 = cmk.Axis(gright[2,1], ylabel="pitch (°)")
ax3 = cmk.Axis(gright[3,1], xlabel="time (s)", ylabel="heading (°)")
cmk.linkxaxes!(ax1, ax2, ax3)


cmk.lines!(ax1, times, roll)
cmk.lines!(ax2, times, pitch)
cmk.lines!(ax3, times, yaw)

f

cmk.colsize!(f.layout, 2, cmk.Auto(0.5))

f



