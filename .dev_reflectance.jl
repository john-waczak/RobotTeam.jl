using RobotTeam
using HDF5
using BenchmarkTools



basepath = "/media/jwaczak/Data/robotteam-data/h5"
test_file = joinpath(basepath, "NoDye_1-1.h5")

@assert isfile(test_file)

generateReflectance!(test_file)


h5 = h5open(test_file, "r")
close(h5)

