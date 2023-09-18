using ..ENVI: HyperspectralImage, read_envi_file

using DataInterpolations
using LoopVectorization
using Images



function generateReflectance!(Reflectance, Radiance, specpath, spechdr, λs)
    # load in downwelling spectrum
    spectrum, hspec, pspec = read_envi_file(specpath, spechdr);
    spec = @view spectrum[:,1,1]

    # load in correction data
    gain, gain_h, gain_ps = read_envi_file(String(gain_path), String(gain_hdr));
    offset, offset_h, offset_ps = read_envi_file(String(offset_path), String(offset_hdr));

    # calculate shutter differences
    cal_shutter = gain_ps["shutter"]
    spec_shutter = pspec["shutter"]
    gain_factor = cal_shutter/spec_shutter # should be a float

    # produce correction frames
    adjusted_gain = gain_factor .* gain[:,1,1]
    frame = spec .- offset[:,1,1]

    # calculate correction
    correction = π .* adjusted_gain .* frame
    clamp!(correction, 0, typemax(eltype(frame)))

    # interpolate the correction to match the datacube's wavelengths
    #----------------------------------------------------------------------
    interp = CubicSpline(correction[:,1,1], pspec["wavelengths"])
    adjustedSpec = interp.(λs)

    # compute reflectance values, clamping them to be ∈ [0, 1]
    @tturbo for λ ∈ axes(Radiance,1), i ∈ axes(Radiance,2), j ∈ axes(Radiance,3)
        Reflectance[λ, i, j] = clamp(π * Radiance[λ, i, j] / adjustedSpec[λ], 0.0, 1.0)
    end
end


function generateReflectance!(Data, specpath, spechdr, λs)
    # load in downwelling spectrum
    spectrum, hspec, pspec = read_envi_file(specpath, spechdr);
    spec = @view spectrum[:,1,1]

    # load in correction data
    gain, gain_h, gain_ps = read_envi_file(String(gain_path), String(gain_hdr));
    offset, offset_h, offset_ps = read_envi_file(String(offset_path), String(offset_hdr));

    # calculate shutter differences
    cal_shutter = gain_ps["shutter"]
    spec_shutter = pspec["shutter"]
    gain_factor = cal_shutter/spec_shutter # should be a float

    # produce correction frames
    adjusted_gain = gain_factor .* gain[:,1,1]
    frame = spec .- offset[:,1,1]

    # calculate correction
    correction = π .* adjusted_gain .* frame
    clamp!(correction, 0, typemax(eltype(frame)))

    # interpolate the correction to match the datacube's wavelengths
    #----------------------------------------------------------------------
    interp = CubicSpline(correction[:,1,1], pspec["wavelengths"])
    adjustedSpec = interp.(λs)

    # compute reflectance values, clamping them to be ∈ [0, 1]
    # @tturbo for j ∈ axes(Data,3), i ∈ axes(Data,2), λ ∈ 1:length(λs)
    #     Data[λ, i, j] = clamp(π * Data[λ, i, j] / adjustedSpec[λ], 0.0, 1.0)
    # end

    Threads.@threads for j ∈ axes(Data,3)
        for i ∈ axes(Data,2)
            for λ ∈ 1:length(λs)
                Data[λ, i, j] = clamp(π * Data[λ, i, j] / adjustedSpec[λ], 0.0, 1.0)
            end
        end
    end

end





