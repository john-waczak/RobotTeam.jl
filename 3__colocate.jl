using RobotTeam

using HDF5
using Dates
using CSV, DataFrames
using HDF5
using ProgressMeter

include("utils/vis_tools.jl")
include("utils/config.jl")

targetpath = "/media/jwaczak/LabData/RobotTeam/prepared/11-23/Targets.csv"
basepath = "/media/jwaczak/LabData/RobotTeam/processed/hsi/11-23"
@assert ispath(targetpath)
@assert ispath(basepath)

df_in = CSV.read(targetpath, DataFrame)
gdf = groupby(df_in, :predye_postdye)


df = gdf[(predye_postdye="Pre-Dye",)]


names(df)

# initialize array and set to NaNs
X_features = Matrix{Float64}(undef, nrow(df), length(features_dict[:varnames]))
X_features .= NaN


function update_features!(X_features, df, collection_id)
    for (root, dirs, files) ∈ walkdir(joinpath(basepath, collection_id))
        for f ∈ files
            if endswith(f, ".h5")
                fpath = joinpath(root, f)
                @info "Working on $(fpath)"
                h5open(fpath, "r") do h5
                    X_h5 = read(h5["data-Δx_0.1/X"])
                    Y_h5 = read(h5["data-Δx_0.1/Y"])
                    IsInbounds = read(h5["data-Δx_0.1/IsInbounds"])

                    Xmin,Xmax = extrema(X_h5)
                    Ymin,Ymax = extrema(Y_h5)

                    Xs = [x for x ∈ X_h5, y ∈ Y_h5]
                    Ys = [y for x ∈ X_h5, y ∈ Y_h5]

                    @showprogress for i ∈ 1:nrow(df)
                        row = @view df[i, :]
                        if row.X ≥ Xmin && row.X ≤ Xmax && row.Y ≥ Ymin && row.Y ≤ Ymax
                            idx = findfirst(row.X .== Xs .&& row.Y .== Ys)
                            if IsInbounds[idx]
                                X_features[i,:] .= h5["data-Δx_0.1/Data_μ"][:,idx]
                            end
                        end
                    end
                end
            end
        end
    end
end



collection_id = "Scotty_2"

update_features!(X_features, df, collection_id)
