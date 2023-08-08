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


# create a plot of the HSI path on a satellite background tile
# w= -97.717472
# n= 33.703572
# s= 33.700797
# e= -97.712413


# satmap = get_background_satmap(w,e,s,n)


# f = cmk.Figure(resolution=(800,600))

# ax = cmk.Axis(f[1,1], xlabel="longitude", ylabel="latitude")

# cmk.heatmap!(ax, satmap.w..satmap.e, satmap.s..satmap.n, satmap.img)
# cmk.xlims!(ax, -97.7168, -97.7125)
# cmk.ylims!(ax, 33.70075, 33.7035)

# f



Δx = 0.1


# reinterpolating the results to a new grid
xmin, xmax = extrema(hsi.X[1,:,:])
ymin, ymax = extrema(hsi.X[2,:,:])


function bump_to_nearest_Δx(val, Δx)
    # units = 0.0:Δx:1.0

    # # assume we are using centimeters

    # whole = floor(val)
    # dec = val-whole

    # dec_new = units[argmin(abs.(units .- dec))]

    # return whole + dec_new
    return round((val ÷ Δx)*Δx, digits=2)
end

bump_to_nearest_Δx(xmin, Δx)


function get_new_bounds(xmin, xmax, ymin, ymax; Δx=0.1)
    Δx_cm = 100*Δx
    @assert 100%Δx_cm == 0

    # pad to nearest Δx
    xmin -= Δx
    xmax += Δx
    ymin -= Δx
    ymax += Δx

   return bump_to_nearest_Δx(xmin, Δx), bump_to_nearest_Δx(xmax, Δx), bump_to_nearest_Δx(ymin, Δx), bump_to_nearest_Δx(ymax, Δx)
end

xmin, xmax, ymin, ymax = get_new_bounds(xmin, xmax, ymin, ymax)


# estimate the bound for minimum pixel spacing in meters
npix = max(size(hsi.X)...)
Δx_min = min((ymax-ymin)/npix, (xmax-xmin)/npix)


# generate a new x-y grid at the desired resolution
xs_new = round.(range(xmin, stop=xmax, step=Δx), digits=2)
ys_new = round.(range(ymin, stop=ymax, step=Δx), digits=2)


Xout = permutedims(cat([x for x ∈ xs_new, y ∈ ys_new], [y for x ∈ xs_new, y ∈ ys_new], dims=3), (3,1,2))  # 1.7 ms
Xhsi = bump_to_nearest_Δx.(hsi.X[1:2,:,:],Δx)

# we should be able to compute the index
findfirst(x->x==Xhsi[1,1,1], xs_new)
findfirst(y->y==Xhsi[2,1,1], ys_new)


# generate hsi pixel indices in outbound grid
Xhsi_is = Matrix{Int}(undef, size(Xhsi, 2), size(Xhsi,3));
Xhsi_js = Matrix{Int}(undef, size(Xhsi, 2), size(Xhsi,3));
@tturbo for j ∈ axes(Xhsi, 3), i ∈ axes(Xhsi,2)
    Xhsi_is[i,j] = Int(round((Xhsi[1,i,j] - xmin)/Δx + 1))
    Xhsi_js[i,j] = Int(round((Xhsi[2,i,j] - ymin)/Δx + 1))
end

# generate boundary mask
IsInbounds = [false for i ∈ axes(Xout,2), j ∈ axes(Xout,3)];
@tturbo for j ∈ axes(Xhsi,3), i ∈ axes(Xhsi,2)
    IsInbounds[Xhsi_is[i,j], Xhsi_js[i,j]] = true
end




# generate array w/ number of HSI pixels per location
Npixels = zeros(Int, size(IsInbounds)...)

@tturbo for j ∈ axes(Xhsi_is,2), i ∈ axes(Xhsi_is, 1)
    Npixels[Xhsi_is[i,j], Xhsi_js[i,j]] += 1
end

cmk.heatmap(Npixels)


# now we should have everything we need to generate the output
