using RobotTeam
using HDF5
import CairoMakie as cmk
using MintsMakieRecipes
cmk.set_theme!(mints_theme)

include("utils/vis_tools.jl")
include("utils/config.jl")


@assert ispath(basepath)
@assert ispath(outpath)
println("Output dir: $outpath")




# create a plot of the HSI path on a satellite background tile
w= -97.717472
n= 33.703572
s= 33.700797
e= -97.712413


satmap = get_background_satmap(w,e,s,n)

for (day, runs) ∈ CollectionsDict
    for (run, fs) ∈ runs
        println("Working on $run")

        savepath = joinpath(outpath, day, run, "maps")
        if !ispath(savepath)
            mkpath(savepath)
        end

        # set up plot
        size_in_inches = (4, 3)
        dpi = 400
        size_in_pixels = size_in_inches .* dpi

        f = cmk.Figure(resolution=size_in_pixels)

        ax = cmk.Axis(f[1,1], xlabel="longitude", ylabel="latitude", title="$day")
        cmk.heatmap!(ax, satmap.w..satmap.e, satmap.s..satmap.n, satmap.img)
        cmk.xlims!(ax, -97.7168, -97.7125)
        cmk.ylims!(ax, 33.70075, 33.7035)

        for f ∈ fs
            fname = split(f.lcfpath, "/")[end-1]
            h5path = joinpath(outpath, day, run, fname * ".h5")
            println("\tplotting $(h5path)")

            try
                h5open(h5path, "r") do h5
                    Img = getRGB(h5; Δx=Δx, α=5, β=0.01)
                    Latitudes = h5["data-Δx_$(Δx)/Latitudes"][:,:]
                    Longitudes = h5["data-Δx_$(Δx)/Longitudes"][:,:]

                    cmk.heatmap!(ax, Longitudes, Latitudes, Img)

                end
            catch e
                println("\n")
                println("FAILED: ", h5path)
                println(e)
                println("\n")
            end
        end

        # save the figure
        save(joinpath(savepath, "map.png"), f)
   end
end



