using RobotTeam
using HDF5
using BenchmarkTools
using DelimitedFiles



# 1. set up test paths

basepath = "/media/jwaczak/Data/robotteam-data/sample-hsi/NoDye_1-1"
outpath = "/media/jwaczak/Data/robotteam-data/h5"
if !ispath(outpath)
    mkpath(outpath)
end


@assert ispath(basepath)

bil_path = joinpath(basepath, "NoDye_1_Pika_XC2_1-radiance.bil")
@assert ispath(bil_path)

bil_hdr = joinpath(basepath, "NoDye_1_Pika_XC2_1-radiance.bil.hdr")
@assert ispath(bil_hdr)

bil_times = joinpath(basepath, "NoDye_1_Pika_XC2_1.bil.times" )
@assert ispath(bil_times)

bil_lcf = joinpath(basepath, "NoDye_1_Pika_XC2_1.lcf")
@assert ispath(bil_lcf)

spec_path = joinpath(basepath, "NoDye_1_downwelling_1_pre.spec")
@assert ispath(spec_path)

spec_hdr = joinpath(basepath, "NoDye_1_downwelling_1_pre.spec.hdr")
@assert ispath(spec_hdr)


# 2. read header

hdr_dict = read_envi_header(bil_hdr)
spec_dict = read_envi_header(spec_hdr)
read_envi_header(bil_times)  # throw error



# 3. get params

get_envi_params(hdr_dict)
get_envi_params(spec_dict)


# 4. read the envi file
@benchmark read_envi_file(bil_path, bil_hdr)     # ~750 ms
@benchmark read_envi_file(spec_path, spec_hdr)   # ~550 μs

envi_to_hdf5(bil_path, bil_hdr, bil_lcf, spec_path, spec_hdr, joinpath(outpath, "NoDye_1-1.h5"))



test_h5 = h5open(joinpath(outpath, "NoDye_1-1.h5"), "r")

read(test_h5["raw/downwelling"], "wavelengths")
read(test_h5["raw/radiance"], "wavelengths")

close(test_h5)
