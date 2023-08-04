module Georectification

using HDF5
using StaticArrays
using LinearAlgebra
using Geodesy
using Meshes
using RelocatableFolders

# set up calibration path
const calibration_path = @path normpath(joinpath(@__DIR__, "../../assets", "calibration"))
const gain_path = @path normpath(joinpath(@__DIR__, "../../assets", "calibration", "gain.spec"))
const gain_hdr = @path normpath(joinpath(@__DIR__, "../../assets", "calibration", "gain.spec.hdr"))
const offset_path = @path normpath(joinpath(@__DIR__, "../../assets", "calibration", "offset.spec"))
const offset_hdr = @path normpath(joinpath(@__DIR__, "../../assets", "calibration", "offset.spec.hdr"))

@assert ispath(calibration_path)
@assert ispath(gain_path)
@assert ispath(gain_hdr)
@assert ispath(offset_path)
@assert ispath(offset_hdr)






include("reflectance.jl")
export generateReflectance!
export generateDerivedMetrics!

include("georectify.jl")
export generateCoords!


include("mesh.jl")
export build_mesh

end
