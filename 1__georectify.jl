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



