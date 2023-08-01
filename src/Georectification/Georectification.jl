module Georectification

using HDF5
using StaticArrays
using LinearAlgebra
using Geodesy


include("reflectance.jl")
export generateReflectance!
export generateDerivedMetrics!

include("georectify.jl")
export generateCoords!



end
