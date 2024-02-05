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


CollectionsDict["12-09"]

bad_hsi_dict = Dict(
    "11-23" => Dict(
        "Scotty_1" => [
            "Scotty_1-16.h6",
            "Scotty_1-26.h5",
            "Scotty_1-27.h5",
            "Scotty_1-28.h5",
            "Scotty_1-29.h5",
            "Scotty_1-30.h5",
            "Scotty_1-31.h5",
            "Scotty_1-32.h5",
        ],
        "Scotty_2" => [
            "Scotty_2-9.h5",
            "Scotty_2-10.h5",
            "Scotty_2-21.h5",
            "Scotty_2-25.h5",
            "Scotty_2-26.h5",
            "Scotty_2-27.h5",
            "Scotty_2-28.h5",
            "Scotty_2-29.h5",
            "Scotty_2-30.h5",
            "Scotty_2-31.h5",
            "Scotty_2-32.h5",
        ],
    ),
    "12-09" => Dict(
        "NoDye_1" => [
            "NoDye_1-1.h5",
            "NoDye_1-9.h5",
            "NoDye_1-19.h5",
            "NoDye_1-23.h5",
        ],
        "NoDye_2" => [
            "NoDye_2-1.h5",
            "NoDye_2-3.h5",
            "NoDye_2-4.h5",
            "NoDye_2-5.h5",
            "NoDye_2-6.h5",
            "NoDye_2-7.h5",
            "NoDye_2-8.h5",
            "NoDye_2-9.h5",
            "NoDye_2-20.h5",
        ],
    ),
    "12-10" => Dict(
        "NoDye_1" => [
            "NoDye_1-14.h5",
            "NoDye_1-19.h5",
            "NoDye_1-20.h5",
            "NoDye_1-23.h5",
        ],
        "NoDye_2" => [
            "NoDye_2-8.h5",
            "NoDye_2-9.h5",
            "NoDye_2-10.h5",
            "NoDye_2-11.h5",
            "NoDye_2-12.h5",
            "NoDye_2-13.h5",
            "NoDye_2-14.h5",
            "NoDye_2-15.h5",
            "NoDye_2-16.h5",
            "NoDye_2-17.h5",
            "NoDye_2-18.h5",
            "NoDye_2-19.h5",
            "NoDye_2-20.h5",
            "NoDye_2-21.h5",
            "NoDye_2-22.h5",
            "NoDye_2-23.h5",
        ]
    ),
)


for (day, runs) ∈ CollectionsDict
    for (run, fs) ∈ runs
        println("Working on $run")

        savepath = joinpath(outpath, day, run, "maps")
        if !ispath(savepath)
            mkpath(savepath)
        end

        # set up plot
        size_in_inches = (8, 3)
        dpi = 400
        size_in_pixels = size_in_inches .* dpi

        fig_main = cmk.Figure(resolution=size_in_pixels)

        ax = cmk.Axis(fig_main[1,1], xlabel="longitude", ylabel="latitude", title="$day")
        ax2 = cmk.Axis(fig_main[1,2])
        cmk.hidespines!(ax2)
        cmk.hidedecorations!(ax2)

        cmk.heatmap!(ax, satmap.w..satmap.e, satmap.s..satmap.n, satmap.img)
        cmk.xlims!(ax, -97.7168, -97.7125)
        cmk.ylims!(ax, 33.70075, 33.7035)

        cmk.xlims!(ax2, -97.7168, -97.7125)
        cmk.ylims!(ax2, 33.70075, 33.7035)

        for f ∈ fs
            fname = split(f.lcfpath, "/")[end-1]
            h5path = joinpath(outpath, day, run, fname * ".h5")
            println("\tplotting $(h5path)")


            fig = cmk.Figure(resolution=(8*400, 3*400))
            ax2_a = cmk.Axis(fig[1,1], xlabel="longitude", ylabel="latitude", title="$day")
            ax2_b = cmk.Axis(fig[1,2])

            cmk.hidespines!(ax2_b)
            cmk.hidedecorations!(ax2_b)

            cmk.heatmap!(ax2_a, satmap.w..satmap.e, satmap.s..satmap.n, satmap.img)
            cmk.xlims!(ax2_a, -97.7168, -97.7125)
            cmk.ylims!(ax2_a, 33.70075, 33.7035)

            try
                h5open(h5path, "r") do h5
                    Img = getRGB(h5; Δx=Δx, α=5, β=0.01)
                    Latitudes = h5["data-Δx_$(Δx)/Latitudes"][:,:]
                    Longitudes = h5["data-Δx_$(Δx)/Longitudes"][:,:]

                    if day in keys(bad_hsi_dict) && run in keys(bad_hsi_dict[day])
                        if !(fname * ".h5" in bad_hsi_dict[day][run])
                            cmk.heatmap!(ax, Longitudes, Latitudes, Img)
                            cmk.heatmap!(ax2, Longitudes, Latitudes, Img)
                        else
                            println("Bad HSI: $(h5path)")
                        end
                    end


                    cmk.heatmap!(ax2_a, Longitudes, Latitudes, Img)
                    cmk.heatmap!(ax2_b, Longitudes, Latitudes, Img)

                    save(joinpath(savepath, fname * ".png"), fig)

                end
            catch e
                println("\n")
                println("FAILED: ", h5path)
                println(e)
                println("\n")
            end
        end

        # save the figure
        save(joinpath(savepath, "map.png"), fig_main)
   end
end

