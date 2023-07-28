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
        ref["wavelengths"] = read(g["radiance/wavelengths"])
    end
end




function generateDerivedMetrics!(h5path::String)
    h5open(h5path, "r+") do h5
        println("generating derived spectral indices")
        g = h5["raw"]
        metrics = create_group(g, "spectral-indices")

        R = read(g["reflectance/reflectance"])
        λs = read(g["reflectance/wavelengths"])

        # mNDWI
        kgreen = findfirst(λs .> 550)
        kswir = findfirst(λs .> 770)
        let
            mNDWI = -1 .* ones(size(R,2), size(R,3))

            green = @view R[kgreen,:,:]
            swir = @view R[kswir,:,:]

            numer = (green .- swir)
            denom = (green .+ swir)

            mNDWI[denom .> 0] = (numer ./ denom)[denom .> 0]

            println("\tsaving mNDWI")
            metrics["mNDWI"] = mNDWI
        end


        # NDVI "normalized difference vegetative index ∈ [-1, 1]"
        kir = findfirst(λs .> 800)
        kred = findfirst(λs .> 680)
        let
            NDVI = -2 .* ones(size(R,2), size(R,3))

            ir_band = @view R[kir,:,:]
            red_band = @view R[kred,:,:]

            numer = (ir_band .- red_band)
            denom = (ir_band .+ red_band)

            NDVI[denom .> 0] .= (numer ./ denom)[denom .> 0]

            println("\tsaving NDVI")
            metrics["NDVI"] = NDVI
        end


        # SR "simple ratio ∈ [0, 30]"
        let
            SR = -1 .* ones(size(R,2), size(R,3))

            ir_band = @view R[kir,:,:]
            red_band = @view R[kred,:,:]
            numer = ir_band
            denom = red_band

            SR[denom .> 0] .= (numer ./ denom)[denom .> 0]

            println("\tsaving SR")
            metrics["SR"] = SR
        end

        # # EVI "enhanced vegetative index ∈ [-1, 1]"
        kblue = findfirst(λs .> 450)
        let
            EVI = -2 .* ones(size(R,2), size(R,3))

            ir_band = @view R[kir,:,:]
            red_band = @view R[kred,:,:]
            blue_band = @view R[kblue,:,:]

            numer = 2.5 .* (ir_band .- red_band)
            denom = ir_band .+ 6 .* red_band .- 7.5 .* blue_band

            EVI[denom .> 0] .= (numer ./ denom)[denom .> 0]

            println("\tsaving EVI")
            metrics["EVI"] = EVI
        end

        # # AVRI "Atmospherical Resistant Vegitative Indes"
        let
            AVRI = -2 .* ones(size(R,2), size(R,3))

            ir_band = @view R[kir,:,:]
            red_band = @view R[kred,:,:]
            blue_band = @view R[kblue,:,:]
            numer = (ir_band .- 2 .* red_band .+ blue_band)
            denom = (ir_band .+ 2 .* red_band .- blue_band)

            AVRI[denom .> 0] .=  (numer ./ denom)[denom .> 0]

            println("\tsaving AVRI")
            metrics["AVRI"] = AVRI
        end


        # # NDVI_705 "Red Edge Normalized Difference Vegetation Index"
        kir = findfirst(λs .> 750)
        kred = findfirst(λs .> 705)
        let
            NDVI_705 = -2 .* ones(size(R,2), size(R,3))

            ir_band = @view R[kir,:,:]
            red_band = @view R[kred,:,:]
            numer = (ir_band .- red_band)
            denom = (ir_band .+ red_band)

            NDVI_705[denom .> 0] .= (numer ./ denom)[denom .> 0]

            println("\tsaving NDVI_705")
            metrics["NDVI_705"] = NDVI_705
        end

        # # MSR_705 "Modified Red Edge Simple Ratio Index"
        kblue = findfirst(λs .> 445)
        let
            MSR_705 = -1 .* ones(size(R,2), size(R,3))

            ir_band = @view R[kir,:,:]
            red_band = @view R[kred,:,:]
            blue_band = @view R[kblue,:,:]

            numer = (ir_band .- blue_band)
            denom = (red_band .- blue_band)

            MSR_705[denom .> 0] .= (numer ./ denom)[denom .> 0]

            println("\tsaving MSR_705")
            metrics["MSR_705"] = MSR_705
        end



        # # MNDVI "modified red edge normalized vegetation index"
        let
            MNDVI = -2 .* ones(size(R,2), size(R,3))

            ir_band = @view R[kir,:,:]
            red_band = @view R[kred,:,:]
            blue_band = @view R[kblue,:,:]

            numer = (ir_band .- red_band)
            denom = (ir_band .+ red_band .- 2 .* blue_band)

            MNDVI[denom .> 0] .=  (numer ./ denom)[denom .> 0]

            println("\tsaving MNDVI")
            metrics["MNDVI"] = MNDVI
        end


        # # VOG1 "vogelmann red edge index"
        kir = findfirst(λs .> 740)
        kred = findfirst(λs .> 720)
        let
            VOG1 = -1 .* ones(size(R,2), size(R,3))

            ir_band = @view R[kir,:,:]
            red_band = @view R[kred,:,:]
            numer = (ir_band)
            denom = (red_band)

            VOG1[denom .> 0] .= (numer ./ denom)[denom .> 0]

            println("\tsaving VOG1")
            metrics["VOG1"] = VOG1
        end


        # # VOG2 "vogelmann red edge index 2"
        k1 = findfirst(λs .> 734)
        k2 = findfirst(λs .> 747)
        k3 = findfirst(λs .> 715)
        k4 = findfirst(λs .> 726)
        let
            VOG2 = -1 .* ones(size(R,2), size(R,3))

            band1 = @view R[k1,:,:]
            band2 = @view R[k2,:,:]
            band3 = @view R[k3,:,:]
            band4 = @view R[k4,:,:]

            numer = (band1 .- band2)
            denom = (band3 .+ band4)

            VOG2[denom .> 0] .=  (numer ./ denom)[denom .> 0]

            println("\tsaving VOG2")
            metrics["VOG2"] = VOG2
        end


        # # VOG3 "vogelmann red edge index 3 ∈ [0, 20]"
        k1 = findfirst(λs .> 734)
        k2 = findfirst(λs .> 747)
        k3 = findfirst(λs .> 715)
        k4 = findfirst(λs .> 720)
        let
            VOG3 = -1 .* ones(size(R,2), size(R,3))

            band1 = @view R[k1,:,:]
            band2 = @view R[k2,:,:]
            band3 = @view R[k3,:,:]
            band4 = @view R[k4,:,:]

            numer = (band1 .- band2)
            denom = (band3 .+ band4)

            VOG3[denom .> 0] .=  (numer ./ denom)[denom .> 0]

            println("\tsaving VOG3")
            metrics["VOG3"] = VOG3
        end


        # # PRI "photochemical reflectance index" ∈[-1, 1]
        k1 = findfirst(λs .> 531)
        k2 = findfirst(λs .> 570)
        let
            PRI = -2 .* ones(size(R,2), size(R,3))

            band1 = @view R[k1,:,:]
            band2 = @view R[k2,:,:]

            numer = (band1 .- band2)
            denom = (band1 .+ band2)

            PRI[denom .> 0] .= (numer ./ denom)[denom .> 0]

            println("\tsaving PRI")

            metrics["PRI"] = PRI
        end

        # # SIPI "structure intensive pigment index" ∈[0, 2]
        k1 = findfirst(λs .> 800)
        k2 = findfirst(λs .> 445)
        k3 = findfirst(λs .> 680)
        let
            SIPI = -1 .* ones(size(R,2), size(R,3))

            band1 = @view R[k1,:,:]
            band2 = @view R[k2,:,:]
            band3 = @view R[k3,:,:]

            numer = (band1 .- band2)
            denom = (band1 .+ band3)

            SIPI[denom .> 0] .= (numer ./ denom)[denom .> 0]

            println("\tsaving SIPI")
            metrics["SIPI"] = SIPI
        end


        # # PSRI "Plant Senescence Reflectance Index"
        k1 = findfirst(λs .> 680)
        k2 = findfirst(λs .> 500)
        k3 = findfirst(λs .> 750)
        let
            PSRI = -2 .* ones(size(R,2), size(R,3))

            band1 = @view R[k1,:,:]
            band2 = @view R[k2,:,:]
            band3 = @view R[k3,:,:]

            numer = (band1 .- band2)
            denom = band3

            PSRI[denom .> 0] .= (numer ./ denom)[denom .> 0]
        end

        # # CRI1 "carotenoid reflectance index"
        k1 = findfirst(λs .> 510)
        k2 = findfirst(λs .> 550)
        let
            CRI1 = -1 .* ones(size(R,2), size(R,3))

            band1 = @view R[k1,:,:]
            band2 = @view R[k2,:,:]

            CRI1[(band1 .> 0) .& (band2 .> 0)] .= ((1 ./ band1) .- (1 ./ band2))[(band1 .> 0) .& (band2 .> 0)]

            println("\tsaving CRI1")
            metrics["CRI1"] = CRI1
        end


        # # CRI2 "carotenoid reflectance index 2"
        k1 = findfirst(λs .> 510)
        k2 = findfirst(λs .> 700)
        let
            CRI2 = -1 .* ones(size(R,2), size(R,3))

            band1 = @view R[k1,:,:]
            band2 = @view R[k2,:,:]

            CRI2[(band1 .> 0) .& (band2 .> 0)] .= ((1 ./ band1) .- (1 ./ band2))[(band1 .> 0) .& (band2 .> 0)]

            println("\tsaving CRI2")
            metrics["CRI2"] = CRI2
        end


        # # ARI1 "anthocyanin reflectance index"
        k1 = findfirst(λs .> 550)
        k2 = findfirst(λs .> 700)
        let
            ARI1 = -1 .* ones(size(R,2), size(R,3))

            band1 = @view R[k1,:,:]
            band2 = @view R[k2,:,:]

            ARI1[(band1 .> 0) .& (band2 .> 0)] = ((1 ./ band1) .- (1 ./ band2))[(band1 .> 0) .& (band2 .> 0)]

            println("\tsaving ARI1")
            metrics["ARI1"] = ARI1
        end


        # # ARI2 "anthocyanin reflectance index 2"
        k1 = findfirst(λs .> 550)
        k2 = findfirst(λs .> 700)
        k3 = findfirst(λs .> 800)
        let
            ARI2 = -1 .* ones(size(R,2), size(R,3))

            band1 = @view R[k1,:,:]
            band2 = @view R[k2,:,:]
            band3 = @view R[k3,:,:]

            ARI2[(band1 .> 0) .& (band2 .> 0)] .= (band3.*((1 ./ band1) .- (1 ./ band2)))[(band1 .> 0) .& (band2 .> 0)]

            println("\tsaving ARI2")
            metrics["ARI2"] = ARI2
        end


        # # WBI "water band index"
        k1 = findfirst(λs .> 900)
        k2 = findfirst(λs .> 970)
        let
            WBI = -1 .* ones(size(R,2), size(R,3))

            band1 = @view R[k1,:,:]
            band2 = @view R[k2,:,:]

            numer = band1
            denom = band2

            WBI[denom .> 0] .= (numer ./ denom)[denom .> 0]

            println("\tsaving WBI")
            metrics["WBI"] = WBI
        end


        # # MCRI "Modified Chlorophyll Absorption Reflectance Index"
        k1 = findfirst(λs .> 550)
        k2 = findfirst(λs .> 670)
        k3 = findfirst(λs .> 701)
        k4 = findfirst(λs .> 780)
        let
            MCRI = -1 .* ones(size(R,2), size(R,3))

            band1 = @view R[k1,:,:]
            band2 = @view R[k2,:,:]
            band3 = @view R[k3,:,:]
            band4 = @view R[k4,:,:]

            MCRI[band2 .> 0] .= (((band3 .-band2) .- 0.2 .* (band3 .- band1)) .* (band3 ./ band2))[band2 .> 0]

            println("\tsaving MCRI")
            metrics["MCRI"] = MCRI
        end


        # # TCARI "transformed chlorophyll absorption reflectance index"
        let
            TCARI = -1 .* ones(size(R,2), size(R,3))

            band1 = @view R[k1,:,:]
            band2 = @view R[k2,:,:]
            band3 = @view R[k3,:,:]

            TCARI[band2 .> 0] .= (3 .* ((band3 .- band2) .- 0.2 .* (band3 .- band1) .* (band3 ./band2)))[band2 .> 0]

            println("\tsaving TCARI")
            metrics["TCARI"] = TCARI
        end
    end
end
