module Georectification

using HDF5
using StaticArrays
using LinearAlgebra

include("reflectance.jl")
export generateReflectance!
export generateDerivedMetrics!

include("georectify.jl")
export generateCoords!



end
