using ..ENVI: read_envi_file
using DataInterpolations
using LoopVectorization



function generateReflectance!(h5path::String; calibration_path::String="calibration")
    gain_path = joinpath(calibration_path, "gain.spec")
    gain_hdr = joinpath(calibration_path, "gain.spec.hdr")
    offset_path = joinpath(calibration_path, "offset.spec")
    offset_hdr = joinpath(calibration_path, "offset.spec.hdr")

    # load in correction data
    gain, gain_h, gain_ps = read_envi_file(gain_path, gain_hdr);
    offset, offset_h, offset_ps = read_envi_file(offset_path, offset_hdr);


    h5open(h5path, "r+") do h5
        # load downwelling and radiance data
        println("\tloading data...")

        spec = read(h5["raw/downwelling"], "irradiance")[:,1,1]
        rad = read(h5["raw/radiance"], "radiance")

        # calculate shutter differences
        println("\tcalculating shutter differences")
        cal_shutter = gain_ps["shutter"]
        spec_shutter = read(h5["raw/downwelling"], "shutter")
        gain_factor = cal_shutter/spec_shutter # should be a float


        # produce correction frames
        println("\tproducing correction frames")
        adjusted_gain = gain_factor .* gain[:,1,1]
        frame = spec .- offset[:,1,1]


        # calculate correction
        println("\tcalculating correction")
        correction = π .* adjusted_gain .* frame
        clamp!(correction, 0, typemax(eltype(frame)))

        # interpolate the correction to match the datacube's wavelengths
        #----------------------------------------------------------------------
        println("\tinterpolating wavelengths")
        interp = CubicSpline(correction[:,1,1], read(h5["raw/downwelling"], "wavelengths"))
        adjustedSpec = interp.(read(h5["raw/radiance"], "wavelengths"))

        # pre-allocate output array
        println("\tpreallocating array")
        R = Array{Float64}(undef, size(h5["raw/radiance/radiance"])...)

        # compute reflectance values, clamping them to be ∈ [0, 1]
        println("\tcomputing reflectances")
        @turbo for λ ∈ axes(rad,1), i ∈ axes(rad,2), j ∈ axes(rad,3)
            R[λ, i, j] = clamp(π * rad[λ, i, j] / adjustedSpec[λ], 0.0, 1.0)
        end


        # write the data to the h5 file
        println("\twriting output to h5 file.")
        g = h5["raw"]
        ref = create_group(g, "reflectance")
        ref["reflectance", chunk=(read(h5["raw/radiance/nbands"]), 1, 1)] = R
    end
end

