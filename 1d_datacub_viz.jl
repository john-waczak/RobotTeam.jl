using RobotTeam
using HDF5
using BenchmarkTools
using CairoMakie
using MintsMakieRecipes

set_theme!(mints_theme)

my_theme = mints_theme
my_theme.Axis.xticklabelsize=20
my_theme.Axis.yticklabelsize=20
my_theme.Axis.xlabelsize=22
my_theme.Axis.ylabelsize=22
my_theme.Axis.titlesize=25
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
λs = read(h5["data-Δx_0.1/λs"])

rgb = getRGB(h5)

close(h5)

indices_metrics = findfirst(x-> x=="mNDWI", varnames):length(varnames)


size_in_inches = (4, 3)
dpi = 300
size_in_pixels = size_in_inches .* dpi

fig = Figure(;resolution=size_in_pixels, figure_padding=50)
ax = CairoMakie.Axis(fig[1,1], xlabel="Longitude", ylabel="Latitude")
hm = heatmap!(ax, satmap.w..satmap.e, satmap.s..satmap.n, satmap.img)

hm_rgb = heatmap!(ax, Longitudes, Latitudes, rgb)

xlims!(ax, -97.717, -97.715)
ylims!(ax, 33.7015, 33.7025)

fig

save("./paper/figures/georectified-hsi-plume.png", fig)

fig


# plot radiance, irradiance, and reflectance for single HSI (12-09 Dye_1-6) for a slice from the plume
fpath = "/media/teamlary/LabData/RobotTeam/raw/hsi/12-09/Dye_1-6"
readdir(fpath)
@assert ispath(fpath)

rad_dat = zeros(UInt16, length(λs))
ref_dat = zeros(Float64, length(λs))
let
    dat, h, p = read_envi_file(joinpath(fpath, "Dye_1_Pika_XC2_6-radiance.bil"), joinpath(fpath, "Dye_1_Pika_XC2_6-radiance.bil.hdr"))
    rad_dat .= dat[:, 100, 100]
end

fnames = joinpath.(fpath, readdir(fpath))

let
    hsi = HyperspectralImage(
        fnames[1],
        fnames[2],
        fnames[4],
        fnames[3],
        fnames[5],
        fnames[6];
        isflipped=true
    );

    ref_dat .= hsi.Reflectance[1:length(λs), 100, 100]
end


irrad_dat = zeros(Float64, 2048)
λs_irrad = zeros(Float64, 2048)
let
    dat, h, p = read_envi_file(fnames[5], fnames[6])
    irrad_dat .= vec(dat)
    λs_irrad .= p["wavelengths"]
end


fig = Figure()
ax = CairoMakie.Axis(fig[1,1], xlabel="λ (nm)", ylabel="Radiance", yticksvisible=false, ygridvisible=false, yminorgridvisible=false, yticklabelsvisible=false)
ax2 = CairoMakie.Axis(fig[1,1], yaxisposition=:right, ylabel="Downwelling Irradiance", yticksvisible=false, ygridvisible=false, yminorgridvisible=false, yticklabelsvisible=false)
linkxaxes!(ax, ax2)

ylims!(ax, 0, nothing)
xlims!(ax, min(λs[1], λs_irrad[1]), max(λs[end], λs_irrad[end]))
ylims!(ax2, 0, nothing)
xlims!(ax2, min(λs[1], λs_irrad[1]), max(λs[end], λs_irrad[end]))


lin = lines!(ax, λs, rad_dat, color=:gray, linewidth=3, label="Radiance")

lin2 = lines!(ax2, λs_irrad, irrad_dat, color=:orange, linewidth=3, label="Downwelling Irradiance")

leg = Legend(fig[1,1], [lin, lin2], ["Radiance", "Downwelling Irradiance"], orientation=:horizontal,  tellheight = false, tellwidth = false, halign=:right, valign=:top, margin = (10, 30, 10, 10),)

fig

save("./paper/figures/radiance-sample.png", fig)
save("./paper/figures/radiance-sample.svg", fig)
save("./paper/figures/radiance-sample.eps", fig)
save("./paper/figures/radiance-sample.pdf", fig)


λ_blue = 495.0
λ_green = 530.0
λ_red = 650.0

cg = cgrad([:purple, :blue, :green, :red, :maroon], [0, λ_blue/1000, λ_green/1000, λ_red/1000  ,1])



fig = Figure()
ax = CairoMakie.Axis(fig[1,1], xlabel="λ (nm)", ylabel="Reflectance", yticksvisible=false, ygridvisible=false, yminorgridvisible=false, yticklabelsvisible=false)

colors = [cg[λ/maximum(λs)] for λ ∈ λs]
lin = lines!(ax, λs, ref_dat, linewidth=3, color=colors)

ylims!(ax, 0, nothing)
xlims!(ax, λs[1], λs[end])

save("./paper/figures/reflectance-sample.png", fig)
save("./paper/figures/reflectance-sample.svg", fig)
save("./paper/figures/reflectance-sample.eps", fig)
save("./paper/figures/reflectance-sample.pdf", fig)
fig


# plot single HSI on map (1209 Dye_1-6)
varnames[indices_metrics[1]]
m1 = Data_μ[indices_metrics[1],:,:]

thresh_h2o = 0.25
idxs_h2o = findall(m1 .≥ thresh_h2o)
idxs_land = findall(m1 .< thresh_h2o)

m1[idxs_land] .= NaN

fig = Figure(;resolution=size_in_pixels);
ax = CairoMakie.Axis(fig[1,1], xlabel="Longitude", ylabel="Latitude", title=printnames[indices_metrics[1]]);
hm = heatmap!(ax, satmap.w..satmap.e, satmap.s..satmap.n, satmap.img)
hm2 = heatmap!(ax, Longitudes, Latitudes, m1, clims=metric_bounds[varnames[indices_metrics[1]]])
fig




# make a map for a single collection 12-09 Dye_1
varnames[indices_metrics[1]]
varnames[indices_metrics[2]]

h2o_thresh = 0.25
idx_mndwi = indices_metrics[1]
idx = indices_metrics[1]
# idx = argmin(abs.(λs .- 630.0))

# set up plot
size_in_inches = (4, 3)
dpi = 400
size_in_pixels = size_in_inches .* dpi
f = Figure(resolution=size_in_pixels)
ax = CairoMakie.Axis(
    f[1, 1],
    xlabel="Longitude",
    ylabel="Latitude",
    xticklabelsize=20,
    yticklabelsize = 20,
    xlabelsize = 25,
    ylabelsize = 25,
    titlesize = 30,
)

heatmap!(ax, satmap.w .. satmap.e, satmap.s .. satmap.n, satmap.img)
xlims!(ax, -97.7168, -97.7125)
ylims!(ax, 33.70075, 33.7035)

f



day = "12-09"
run = "Dye_1"

for f ∈ CollectionsDict[day][run]
    fname = split(f.lcfpath, "/")[end-1]

    if !(fname ∈ [
        run * "-20",
        run * "-21",
        run * "-22",
    ])

        h5path = joinpath(outpath, day, run, fname * ".h5")
        println("\t\tplotting $(h5path)")
        try
            h5open(h5path, "r") do h5
                Latitudes = h5["data-Δx_$(Δx)/Latitudes"][:, :]
                Longitudes = h5["data-Δx_$(Δx)/Longitudes"][:, :]
                Data_mndwi = h5["data-Δx_$(Δx)/Data_μ"][idx_mndwi, :, :]
                Data = h5["data-Δx_$(Δx)/Data_μ"][idx, :, :]

                # get indices of non-water pixels
                idx_not_h2o = findall(Data_mndwi .< h2o_thresh)
                # set value to NaN for non-water
                Data[idx_not_h2o] .= NaN

                heatmap!(ax, Longitudes, Latitudes, Data; clims=metric_bounds[varnames[idx]])
                # heatmap!(ax, Longitudes, Latitudes, Data; clims=(0, 1))
                f
            end
        catch e
          println("\n")
          println("FAILED: ", h5path)
          println(e)
          println("\n")
        end
    end
end


# create colorbar
cb = Colorbar(f[1,2], limits=metric_bounds[varnames[idx]], label="$(printnames[idx])", labelsize=25, ticklabelsize=20)

f

save("paper/figures/mNDWI_example.png", f)



# loop over all data and create maps

# # loop through and make maps
# for (day, runs) ∈ CollectionsDict
#     for (run, fs) ∈ runs
#         println("Working on $run")

#         savepath = joinpath(outpath, day, run, "maps")
#         if !ispath(savepath)
#             mkpath(savepath)
#         end

#         # loop over all metrics
#         for idx ∈ indices_metrics
#             GC.gc()

#             println("\t$(printnames[idx])")

#             # set up plot
#             size_in_inches = (4, 3)
#             dpi = 400
#             size_in_pixels = size_in_inches .* dpi
#             f = Figure(resolution=size_in_pixels)
#             ax = CairoMakie.Axis(f[1,1], xlabel="longitude", ylabel="latitude", title="$day")
#             heatmap!(ax, satmap.w..satmap.e, satmap.s..satmap.n, satmap.img)
#             xlims!(ax, -97.7168, -97.7125)
#             ylims!(ax, 33.70075, 33.7035)

#             # loop over each file in the run and add the values
#             for f ∈ fs
#                 fname = split(f.lcfpath, "/")[end-1]
#                 h5path = joinpath(outpath, day, run, fname * ".h5")
#                 println("\t\tplotting $(h5path)")

#                 try
#                     h5open(h5path, "r") do h5
#                         Latitudes = h5["data-Δx_$(Δx)/Latitudes"][:,:]
#                         Longitudes = h5["data-Δx_$(Δx)/Longitudes"][:,:]
#                         Data = h5["data-Δx_$(Δx)/Data_μ"][idx,:,:]

#                         heatmap!(ax, Longitudes, Latitudes, Data; clims=metric_bounds[varnames[idx]])
#                         f
#                     end
#                 catch e
#                     println("\n")
#                     println("FAILED: ", h5path)
#                     println(e)
#                     println("\n")
#                 end
#             end

#             cb = Colorbar(f[1,2], limits=metric_bounds[varnames[idx]], label="$(printnames[idx])")

#             # save the figure
#             png_path = joinpath(savepath, "map_$(varnames[idx]).png")
#             save(png_path, f)
#         end

#    end
# end



