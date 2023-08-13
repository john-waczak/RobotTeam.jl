basepath = "/media/jwaczak/LabData/RobotTeam/raw/hsi"
@assert ispath(basepath)

outpath = "/media/jwaczak/LabData/RobotTeam/processed/hsi"
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


# define bounding box for Scotty's Ranch
w= -97.717472
n= 33.703572
s= 33.700797
e= -97.712413

