using RobotTeam

using Dates
using CSV, DataFrames
using HDF5

include("utils/vis_tools.jl")
include("utils/config.jl")

boatpath = "/media/teamlary/LabData/RobotTeam/raw/boat"
@assert ispath(boatpath)

outpath = "/media/teamlary/LabData/RobotTeam/processed/boat"
if !ispath(outpath)
    mkpath(outpath)
end

collection_dates = [
    ["20201123", "11-23"],
    ["20201209", "12-09"],
    ["20201210", "12-10"],
    ["20220324", "03-24"],
]


paths_dict = Dict()
for ds ∈ collection_dates
    din, dout = ds
    paths_dict[dout] = (;raw=joinpath(boatpath,din), processed=joinpath(outpath,dout))
    if !ispath(joinpath(outpath, dout))
        mkpath(joinpath(outpath, dout))
    end
end



for (date, paths) ∈ paths_dict
    @info "Processing $(paths.raw)"
    processBoatFiles(paths.raw, paths.processed)
end


prepared_path = "/media/teamlary/LabData/RobotTeam/prepared"
if !ispath(prepared_path)
    mkpath(prepared_path)
    for dts ∈ collection_dates
        mkpath(joinpath(prepared_path, dts[2]))
    end
end



dirpath = paths_dict["11-23"].processed
prepared_path

for (day, fs) ∈ paths_dict
    dirpath = fs.processed
    println("Working on $dirpath")
    try
        combine_boat_dfs(dirpath, joinpath(prepared_path, day); Δx=Δx)
    catch e
        println("FAILED: ", dirpath)
        println(e)
    end
end


# now we want to get the "category" for the flight so we can carefully sort by time.
hsi_path = "/media/teamlary/LabData/RobotTeam/raw/hsi"
@assert ispath(hsi_path)



collection_times = Dict()

for (day, collections) ∈ CollectionsDict
    if day != "03-24"
        collection_times[day] = Dict()
        for collection ∈ keys(collections)
            tstarts = []
            tends = []
            hsi_fs = get_raw_file_list.(get_bil_files(joinpath(hsi_path, day), collection))
            for hsi_f ∈ hsi_fs
                fd = FlightData(hsi_f.lcfpath, hsi_f.timespath)
                
                push!(tstarts, fd.start_time)
                push!(tends, fd.start_time + Millisecond(round(Int, fd.times[end]*1000 + 1)))
            end

            collection_times[day][collection] = (; tstart=minimum(tstarts), tend=maximum(tends))
        end
    end
end

collection_times["11-23"]["Scotty_1"]

predye_postdye_dict = Dict(
    "11-23" => (;
                tbefore = DateTime(2020, 11, 23, 17, 56),
                t0 = DateTime(2020, 11, 23, 18, 53),
                t1 = DateTime(2020, 11, 23, 19, 20)),
    "12-09" => (;
                tbefore = DateTime(2020, 12, 09, 13, 47, 44),
                t0 = DateTime(2020, 12, 09, 14, 40, 00),
                t1 = DateTime(2020, 12, 09, 14, 54, 00),
                ),
    "12-10" => (;
                tbefore = DateTime(2020, 12, 10, 15, 18, 00),
                t0 = DateTime(2020, 12, 10, 18, 13, 00),
                t1 = DateTime(2020, 12, 10, 18, 43, 00),
                )

)


for (day, categories) ∈ collection_times
    df_path = joinpath(prepared_path, day, "Targets.csv")
    @assert ispath(df_path)
    df_targets = CSV.File(df_path) |> DataFrame
    df_targets.category = ["PostFlights" for _ ∈ 1:nrow(df_targets)]
    df_targets.predye_postdye = ["ignore" for _ ∈ 1:nrow(df_targets)]

    cats = [key for key ∈ keys(collection_times[day])]
    tstarts = [collection_times[day][cat].tstart for cat ∈ cats]
    tends = [collection_times[day][cat].tend for cat ∈ cats]
    df_cats = DataFrame(cats=cats, tstarts=tstarts, tends=tends)
    sort!(df_cats, :tstarts)

    for i ∈ 1:nrow(df_targets)
        t = df_targets.utc_dt[i]


        # add the flight category
        for j ∈ nrow(df_cats):-1:1
            if t ≤ df_cats.tends[j]
                df_targets.category[i] = df_cats.cats[j]
            end

            if t < df_cats.tstarts[j]
                df_targets.category[i] = df_cats.cats[j]*"_preflight"
            end
        end

        # add the category for pre/post dye release
        if t ≥ predye_postdye_dict[day].tbefore && t ≤ predye_postdye_dict[day].t0
            df_targets.predye_postdye[i] = "Pre-Dye"
        elseif t ≥ predye_postdye_dict[day].t0 && t ≤ predye_postdye_dict[day].t1
            df_targets.predye_postdye[i] = "Post_Dye"
        end

    end

    # df_targets.category
    # df_targets.predye_postdye
    CSV.write(joinpath(prepared_path, day, "Targets.csv"), df_targets)
end


