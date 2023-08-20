using RobotTeam

using HDF5
using Dates
using CSV, DataFrames
using HDF5
using ProgressMeter

include("utils/vis_tools.jl")
include("utils/config.jl")

outpath = "/media/teamlary/LabData/RobotTeam/finalized"
if !ispath(outpath)
    mkpath(outpath)
    mkpath(joinpath(outpath, "11-23"))
    mkpath(joinpath(outpath, "12-09"))
    mkpath(joinpath(outpath, "12-10"))
end


function update_features!(X_features, df, collection_id)
    for (root, dirs, files) ∈ walkdir(joinpath(basepath, collection_id))
        for f ∈ files
            if endswith(f, ".h5")
                fpath = joinpath(root, f)
                println("Working on $(fpath)")
                h5open(fpath, "r") do h5
                    X_h5 = read(h5["data-Δx_0.1/X"])
                    Y_h5 = read(h5["data-Δx_0.1/Y"])
                    IsInbounds = read(h5["data-Δx_0.1/IsInbounds"])

                    Xmin,Xmax = extrema(X_h5)
                    Ymin,Ymax = extrema(Y_h5)

                    Data = read(h5["data-Δx_0.1/Data_μ"])

                    Xs = [x for x ∈ X_h5, y ∈ Y_h5]
                    Ys = [y for x ∈ X_h5, y ∈ Y_h5]

                    for i ∈ 1:nrow(df)
                        if all(isnan.(X_features[i,:]))
                            row = @view df[i, :]
                            if row.X ≥ Xmin && row.X ≤ Xmax && row.Y ≥ Ymin && row.Y ≤ Ymax
                                idx = findfirst(row.X .== Xs .&& row.Y .== Ys)
                                if IsInbounds[idx]
                                    X_features[i,:] .= Data[:,idx]
                                end
                            end
                        end
                    end

                end
            end
        end
    end
end



# now lets go through in order:
# 11-23

targetpath = "/media/teamlary/LabData/RobotTeam/prepared/11-23/Targets.csv"
basepath = "/media/teamlary/LabData/RobotTeam/processed/hsi/11-23"
@assert ispath(targetpath)
@assert ispath(basepath)

df_in = CSV.read(targetpath, DataFrame)
gdf = groupby(df_in, :predye_postdye)

df = gdf[(predye_postdye="Pre-Dye",)]

names(df)

# initialize array and set to NaNs
X_features = Matrix{Float64}(undef, nrow(df), length(features_dict[:varnames]))
X_features .= NaN


# first let's start with collection 2 since it's closer to the boat time
update_features!(X_features, df, "Scotty_2")

idxs_nans = [all(isnan.(X_features[i,:])) for i in 1:size(X_features, 1)]
println(sum(idxs_nans))

update_features!(X_features, df, "Scotty_1")

idxs_nans = [all(isnan.(X_features[i,:])) for i in 1:size(X_features, 1)]
println(sum(idxs_nans))

idxs_notnans = [!all(isnan.(X_features[i,:])) for i in 1:size(X_features, 1)]

df_features = DataFrame(X_features[idxs_notnans, :], features_dict[:varnames])
df_targets = df[idxs_notnans,:]

@assert nrow(df_targets) == nrow(df_features)
# save results:
CSV.write(joinpath(outpath, "11-23", "Targets.csv"), df_targets)
CSV.write(joinpath(outpath, "11-23", "Features.csv"), df_features)



# repeat for 12-09
targetpath = "/media/teamlary/LabData/RobotTeam/prepared/12-09/Targets.csv"
basepath = "/media/teamlary/LabData/RobotTeam/processed/hsi/12-09"
@assert ispath(targetpath)
@assert ispath(basepath)

df_in = CSV.read(targetpath, DataFrame)
gdf = groupby(df_in, :predye_postdye)

df = gdf[(predye_postdye="Pre-Dye",)]

names(df)
df.category

# initialize array and set to NaNs
X_features = Matrix{Float64}(undef, nrow(df), length(features_dict[:varnames]))
X_features .= NaN


update_features!(X_features, df, "NoDye_1")

idxs_nans = [all(isnan.(X_features[i,:])) for i in 1:size(X_features, 1)]
println(sum(idxs_nans))

update_features!(X_features, df, "NoDye_2")

idxs_nans = [all(isnan.(X_features[i,:])) for i in 1:size(X_features, 1)]
println(sum(idxs_nans))

idxs_notnans = [!all(isnan.(X_features[i,:])) for i in 1:size(X_features, 1)]

df_features = DataFrame(X_features[idxs_notnans, :], features_dict[:varnames])
df_targets = df[idxs_notnans,:]

@assert nrow(df_targets) == nrow(df_features)
# save results:
CSV.write(joinpath(outpath, "12-09", "Targets.csv"), df_targets)
CSV.write(joinpath(outpath, "12-09", "Features.csv"), df_features)


# repeat for 12-10
targetpath = "/media/teamlary/LabData/RobotTeam/prepared/12-10/Targets.csv"
basepath = "/media/teamlary/LabData/RobotTeam/processed/hsi/12-10"
@assert ispath(targetpath)
@assert ispath(basepath)

df_in = CSV.read(targetpath, DataFrame)
gdf = groupby(df_in, :predye_postdye)

df = gdf[(predye_postdye="Pre-Dye",)]

names(df)
df.category

# initialize array and set to NaNs
X_features = Matrix{Float64}(undef, nrow(df), length(features_dict[:varnames]))
X_features .= NaN


update_features!(X_features, df, "NoDye_1")

idxs_nans = [all(isnan.(X_features[i,:])) for i in 1:size(X_features, 1)]
println(sum(idxs_nans))

update_features!(X_features, df, "NoDye_2")

idxs_nans = [all(isnan.(X_features[i,:])) for i in 1:size(X_features, 1)]
println(sum(idxs_nans))

idxs_notnans = [!all(isnan.(X_features[i,:])) for i in 1:size(X_features, 1)]

df_features = DataFrame(X_features[idxs_notnans, :], features_dict[:varnames])
df_targets = df[idxs_notnans,:]

@assert nrow(df_targets) == nrow(df_features)
# save results:
CSV.write(joinpath(outpath, "12-10", "Targets.csv"), df_targets)
CSV.write(joinpath(outpath, "12-10", "Features.csv"), df_features)


