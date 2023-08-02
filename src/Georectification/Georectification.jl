module Georectification

using HDF5
using StaticArrays
using LinearAlgebra
using Geodesy
using Meshes

include("reflectance.jl")
export generateReflectance!
export generateDerivedMetrics!

include("georectify.jl")
export generateCoords!


include("mesh.jl")
export build_mesh

end
