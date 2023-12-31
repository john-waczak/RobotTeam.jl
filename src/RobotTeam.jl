module RobotTeam


include("ENVI/ENVI.jl")
# ENVI.jl exports
using .ENVI
export FileNotAnEnviHeader
export EnviHeaderParsingError
export read_envi_header
export get_envi_params
export read_envi_file
export FlightData, HyperspectralImage


include("Georectification/Georectification.jl")
using .Georectification

export generateReflectance!, generateReflectance2!
export generateDerivedMetrics!
export generateCoords!
export bump_to_nearest_Δx, get_new_bounds, get_resampled_grid, resample_datacube, save_resampled_hsi
export resample_datacube_fast, save_resampled_hsi_fast
export generate_derived_metrics!
export get_bil_files, get_raw_file_list



include("Boat/Boat.jl")
using .Boat
export importAirMar, importCOM1, importCOM2, importCOM3, importLISST, importNMEA, processBoatFiles, combine_boat_dfs

# visualization tools

include("visualization.jl")
using .Visualization
export process_image, getRGB, vis_rectified_cube

end
