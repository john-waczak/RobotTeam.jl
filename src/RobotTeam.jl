module RobotTeam

include("ENVI/ENVI.jl")

# ENVI.jl exports
using .ENVI
export FileNotAnEnviHeader
export EnviHeaderParsingError
export read_envi_header
export get_envi_params
export read_envi_file
export envi_to_hdf5

include("Georectification/Georectification.jl")
using .Georectification

export generateReflectance!
export generateDerivedMetrics!
export generateCoords!


end
