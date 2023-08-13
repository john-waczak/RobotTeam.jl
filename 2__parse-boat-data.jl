using RobotTeam

using CSV, DataFrames

include("utils/vis_tools.jl")
include("utils/config.jl")

boatpath = "/media/jwaczak/LabData/RobotTeam/raw/boat"
@assert ispath(boatpath)

outpath = "/media/jwaczak/LabData/RobotTeam/processed/boat"
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


prepared_path = "/media/jwaczak/LabData/RobotTeam/prepared"
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
