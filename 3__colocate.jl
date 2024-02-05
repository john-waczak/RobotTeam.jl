using RobotTeam

using Statistics
using HDF5
using Dates
using CSV, DataFrames
using HDF5
using ProgressMeter

include("utils/vis_tools.jl")
include("utils/config.jl")

outpath = "/Users/johnwaczak/data/robot-team/finalized"
if !ispath(outpath)
    mkpath(outpath)
    mkpath(joinpath(outpath, "11-23"))
    mkpath(joinpath(outpath, "12-09"))
    mkpath(joinpath(outpath, "12-10"))
    mkpath(joinpath(outpath, "Full"))
end



function sort_files_numerical!(flist)
    idx_sort = sortperm([lpad(split(s[1], "-")[end], 2, "0") for s in split.(flist, ".")])
    flist .= flist[idx_sort]
end


basepath_11_23 = "/Users/johnwaczak/data/robot-team/processed/hsi/11-23"
basepath_12_09 = "/Users/johnwaczak/data/robot-team/processed/hsi/12-09"
basepath_12_10 = "/Users/johnwaczak/data/robot-team/processed/hsi/12-10"


df_11_23 = CSV.read("/Users/johnwaczak/data/robot-team/prepared/11-23/Targets.csv", DataFrame);
df_12_09 = CSV.read("/Users/johnwaczak/data/robot-team/prepared/12-09/Targets.csv", DataFrame);
df_12_10 = CSV.read("/Users/johnwaczak/data/robot-team/prepared/12-10/Targets.csv", DataFrame);


collections = Dict(
    "11-23" => Dict(
        "Scotty_1" => sort_files_numerical!([joinpath(basepath_11_23, "Scotty_1", f) for f in readdir(joinpath(basepath_11_23, "Scotty_1")) if endswith(f, ".h5")]),
        "Scotty_2" => sort_files_numerical!([joinpath(basepath_11_23, "Scotty_2", f) for f in readdir(joinpath(basepath_11_23, "Scotty_2")) if endswith(f, ".h5")])
    ),
    "12-09" => Dict(
        "NoDye_1" => sort_files_numerical!([joinpath(basepath_12_09, "NoDye_1", f) for f in readdir(joinpath(basepath_12_09, "NoDye_1")) if endswith(f, ".h5")]),
        "NoDye_2" => sort_files_numerical!([joinpath(basepath_12_09, "NoDye_2", f) for f in readdir(joinpath(basepath_12_09, "NoDye_2")) if endswith(f, ".h5")]),
    ),
    "12-10" => Dict(
        "NoDye_1" => sort_files_numerical!([joinpath(basepath_12_10, "NoDye_1", f) for f in readdir(joinpath(basepath_12_10, "NoDye_1")) if endswith(f, ".h5")]),
        "NoDye_2" => sort_files_numerical!([joinpath(basepath_12_10, "NoDye_2", f) for f in readdir(joinpath(basepath_12_10, "NoDye_2")) if endswith(f, ".h5")]),
    )
)


function find_matching_data(h5_path, df)
    df_features_out = DataFrame()
    df_targets_out = DataFrame()

    h5open(h5_path, "r") do h5
        # read position in meters
        X_h5 = read(h5["data-Δx_$(Δx)/X"])
        Y_h5 = read(h5["data-Δx_$(Δx)/Y"])
        varnames_h5 = read(h5["data-Δx_$(Δx)/varnames"])
        @assert all(features_dict[:varnames] .== varnames_h5)
        IsInbounds = read(h5["data-Δx_$(Δx)/IsInbounds"])

        # get bounds of datacube
        Xmin,Xmax = extrema(X_h5)
        Ymin,Ymax = extrema(Y_h5)

        # read in the data
        Data = read(h5["data-Δx_$(Δx)/Data"])

        # generate coordinate matrices of same dim as Data
        Xs = [x for x ∈ X_h5, y ∈ Y_h5]
        Ys = [y for x ∈ X_h5, y ∈ Y_h5]


        # loop over the rows of our boat data to find all matching boat data
        idx_match_hsi = CartesianIndex{2}[]
        idx_match_df = Int[]

        for i ∈ 1:nrow(df)
            row = @view df[i, :]
            # filter by bounding box
            if row.X ≥ Xmin && row.X ≤ Xmax && row.Y ≥ Ymin && row.Y ≤ Ymax
                idx = findfirst(row.X .== Xs .&& row.Y .== Ys)
                if !isnothing(idx)

                    # once we have a match, create a mask surrouding
                    # the pixel to obtain a 3×3 grid to account
                    # for the size of the boat and length
                    # of the sampling tubes...
                    idxs = [idx - CartesianIndex(ii,jj) for ii in -1:1 for jj in -1:1 if checkbounds(Bool, Xs, (idx - CartesianIndex(ii,jj)))]

                    for idx_good in idxs
                        if IsInbounds[idx_good]
                            push!(idx_match_hsi, idx_good)
                            push!(idx_match_df, i)
                        end
                    end

                    # if IsInbounds[idx]
                    #     # push!(df_features_out, Data[:,idx])
                    #     # push!(df_targets_out, df[i,:])
                    #     push!(idx_match_hsi, idx)
                    #     push!(idx_match_df, i)
                    # end
                end
            end
        end

        df_features_out = DataFrame(Data[:,idx_match_hsi]', features_dict[:varnames])
        df_targets_out = df[idx_match_df,:]
    end

    return df_features_out, df_targets_out
end


function get_all_matching_data(collection, collection_id)
    df_in = CSV.read(joinpath("/Users/johnwaczak/data/robot-team/prepared/$(collection)/Targets.csv"), DataFrame)
    gdf = groupby(df_in, :predye_postdye)
    df = gdf[(predye_postdye="Pre-Dye",)]

    if collection == "12-10"
        df = df[df.category .== "Dye_1_preflight", :]  # add this to remove the other noisy points that aren't helping
    end


    # data_features = DataFrame[]
    # data_targets = DataFrame[]

    @showprogress for h5_path in collections[collection][collection_id]
        fname = split(split(h5_path, "/")[end], ".")[1]

        df_f, df_t = find_matching_data(h5_path, df)

        savepath = joinpath(outpath, collection, collection_id)
        if !ispath(savepath)
            mkpath(savepath)
        end

        CSV.write(joinpath(savepath, fname*"__Features.csv"), df_f)
        CSV.write(joinpath(savepath, fname*"__Targets.csv"), df_t)

        # push!(data_features, df_f)
        # push!(data_targets, df_t)
    end

    # return vcat(data_features...), vcat(data_targets...)
end



# use CDOM as a proxy to filter values, keeping 95% of data
function idx_filter(df_features, df_targets)
    q_low = quantile(df_targets.CDOM, 0.025)
    q_high = quantile(df_targets.CDOM, 0.975)
    idx_filtered = findall(df_targets.CDOM .>= q_low .&& df_targets.CDOM .<= q_high)
    return idx_filtered
end



get_all_matching_data("11-23", "Scotty_1")
get_all_matching_data("11-23", "Scotty_2")
get_all_matching_data("12-09", "NoDye_1")
get_all_matching_data("12-09", "NoDye_2")
get_all_matching_data("12-10", "NoDye_1")
get_all_matching_data("12-10", "NoDye_2")





# # process 11-23
# df_features_11_23_1, df_targets_11_23_1 = get_all_matching_data("11-23", "Scotty_1")
# # filter to sane values
# idx_filtered = idx_filter(df_features_11_23_1, df_targets_11_23_1)
# df_features_11_23_1 = df_features_11_23_1[idx_filtered,:]
# df_targets_11_23_1 = df_targets_11_23_1[idx_filtered,:]

# @assert nrow(df_features_11_23_1) == nrow(df_targets_11_23_1)
# CSV.write(joinpath(outpath, "11-23", "Targets_1.csv"), df_targets_11_23_1)
# CSV.write(joinpath(outpath, "11-23", "Features_1.csv"), df_features_11_23_1)


# df_features_11_23_2, df_targets_11_23_2 = get_all_matching_data("11-23", "Scotty_2")
# # filter to sane values
# idx_filtered = idx_filter(df_features_11_23_2, df_targets_11_23_2)
# df_features_11_23_2 = df_features_11_23_2[idx_filtered,:]
# df_targets_11_23_2 = df_targets_11_23_2[idx_filtered,:]

# @assert nrow(df_features_11_23_1) == nrow(df_targets_11_23_1)
# CSV.write(joinpath(outpath, "11-23", "Targets_2.csv"), df_targets_11_23_2)
# CSV.write(joinpath(outpath, "11-23", "Features_2.csv"), df_features_11_23_2)


# # process 12-09
# df_features_12_09_1, df_targets_12_09_1 = get_all_matching_data("12-09", "NoDye_1")
# # filter to sane values
# idx_filtered = idx_filter(df_features_12_09_1, df_targets_12_09_1)
# df_features_12_09_1 = df_features_12_09_1[idx_filtered,:]
# df_targets_12_09_1 = df_targets_12_09_1[idx_filtered,:]

# @assert nrow(df_features_12_09_1) == nrow(df_targets_12_09_1)
# CSV.write(joinpath(outpath, "12-09", "Targets_1.csv"), df_targets_12_09_1)
# CSV.write(joinpath(outpath, "12-09", "Features_1.csv"), df_features_12_09_1)

# df_features_12_09_2, df_targets_12_09_2 = get_all_matching_data("12-09", "NoDye_2")
# # filter to sane values
# idx_filtered = idx_filter(df_features_12_09_2, df_targets_12_09_2)
# df_features_12_09_2 = df_features_12_09_2[idx_filtered,:]
# df_targets_12_09_2 = df_targets_12_09_2[idx_filtered,:]

# @assert nrow(df_features_12_09_2) == nrow(df_targets_12_09_2)
# CSV.write(joinpath(outpath, "12-09", "Targets_2.csv"), df_targets_12_09_2)
# CSV.write(joinpath(outpath, "12-09", "Features_2.csv"), df_features_12_09_2)


# # process 12-10
# df_features_12_10_1, df_targets_12_10_1 = get_all_matching_data("12-10", "NoDye_1")
# # filter to sane values
# idx_filtered = idx_filter(df_features_12_10_1, df_targets_12_10_1)
# df_features_12_10_1 = df_features_12_10_1[idx_filtered,:]
# df_targets_12_10_1 = df_targets_12_10_1[idx_filtered,:]

# @assert nrow(df_features_12_10_1) == nrow(df_targets_12_10_1)
# CSV.write(joinpath(outpath, "12-10", "Targets_1.csv"), df_targets_12_10_1)
# CSV.write(joinpath(outpath, "12-10", "Features_1.csv"), df_features_12_10_1)

# df_features_12_10_2, df_targets_12_10_2 = get_all_matching_data("12-10", "NoDye_2")
# # filter to sane values
# idx_filtered = idx_filter(df_features_12_10_2, df_targets_12_10_2)
# df_features_12_10_2 = df_features_12_10_2[idx_filtered,:]
# df_targets_12_10_2 = df_targets_12_10_2[idx_filtered,:]

# @assert nrow(df_features_12_10_2) == nrow(df_targets_12_10_2)
# CSV.write(joinpath(outpath, "12-10", "Targets_2.csv"), df_targets_12_10_2)
# CSV.write(joinpath(outpath, "12-10", "Features_2.csv"), df_features_12_10_2)


# # join together for full dataset
# df_features_full_1 = vcat(
#     df_features_11_23_1,
#     df_features_12_09_1,
#     df_features_12_10_1,
# )

# df_targets_full_1 = vcat(
#     df_targets_11_23_1,
#     df_targets_12_09_1,
#     df_targets_12_10_1,
# )

# @assert nrow(df_features_full_1) == nrow(df_targets_full_1)
# CSV.write(joinpath(outpath, "Full", "Targets_1.csv"), df_targets_full_1)
# CSV.write(joinpath(outpath, "Full", "Features_1.csv"), df_features_full_1)



# df_features_full_2 = vcat(
#     df_features_11_23_2,
#     df_features_12_09_2,
#     df_features_12_10_2,
# )

# df_targets_full_2 = vcat(
#     df_targets_11_23_2,
#     df_targets_12_09_2,
#     df_targets_12_10_2,
# )

# @assert nrow(df_features_full_2) == nrow(df_targets_full_2)
# CSV.write(joinpath(outpath, "Full", "Targets_2.csv"), df_targets_full_2)
# CSV.write(joinpath(outpath, "Full", "Features_2.csv"), df_features_full_2)


