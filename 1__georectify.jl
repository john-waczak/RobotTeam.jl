using RobotTeam
using HDF5
import CairoMakie as cmk
using MintsMakieRecipes
cmk.set_theme!(mints_theme)

include("utils/vis_tools.jl")


basepath = "/media/jwaczak/LabData/RobotTeam/raw/hsi"
@assert ispath(basepath)

outpath = "/media/jwaczak/LabData/RobotTeam/processed"
if !ispath(outpath)
    mkpath(outpath)
end


CollectionsDict = Dict(
    "11-23" => Dict(
        "Scotty_1" => [],
        "Scotty_2" => [],
        "Scotty_3" => [],
        "Scotty_4" => [],
        "Scotty_5" => [],
    ),
    "12-09" => Dict(
        "NoDye_1" => [],
        "NoDye_2" => [],
        "Dye_1" => [],
        "Dye_2" => [],
    ),
    "12-10" => Dict(
        "NoDye_1" => [],
        "NoDye_2" => [],
        "Dye_1" => [],
        "Dye_2" => [],
    ),
    "03-24" => Dict(
        "Demonstration" => [],
        "Demonstration_long" => [],
    )
)


# look up bil files
for (day, runs) ∈ CollectionsDict
    for run ∈ keys(runs)
        if !ispath(joinpath(outpath, day, run))
            mkpath(joinpath(outpath, day, run))
        end

        CollectionsDict[day][run] = get_raw_file_list.(get_bil_files(joinpath(basepath, day), run))
    end
end







θ_view=30.8              # viewing angle
z_ground=292.0           # ground height (m)
isflipped=false          # flip pixel orientation
Δx = 0.10                # resampling resolution (m)
is_spec_chunked=false    # chunk HDF5 by pixel
is_band_chunked=true     # chunk HDF5 by band


CollectionsDict

for (day, runs) ∈ CollectionsDict
    for (run, fs) ∈ runs
        for f ∈ fs
            println("Processing $(f.bilpath)")
            try

                println("\topening HSI")
                hsi = HyperspectralImage(
                    f.bilpath,
                    f.bilhdr,
                    f.lcfpath,
                    f.timespath,
                    f.specpath,
                    f.spechdr;
                    isflipped=true
                )

                println("\tresampling to new grid")
                xs, ys, isnorth, zone, Longitudes, Latitudes, IsInbounds, varnames, printnames, λs, Data_μ, Data_σ = resample_datacube(hsi; Δx=Δx)


                fname = split(f.lcfpath, "/")[end-1] * ".h5"

                println("\tsaving as $fname")
                save_resampled_hsi(
                    xs,
                    ys,
                    isnorth,
                    zone,
                    Longitudes,
                    Latitudes,
                    IsInbounds,
                    varnames,
                    printnames,
                    λs,
                    Data_μ,
                    Data_σ,
                    joinpath(outpath, day, run, fname);
                    Δx=Δx,
                    is_spec_chunked=is_spec_chunked,
                    is_band_chunked=is_band_chunked
                )

            catch e
                println("FAILED: $(f.bilpath)")
                println(e)
            end
        end
    end
end



