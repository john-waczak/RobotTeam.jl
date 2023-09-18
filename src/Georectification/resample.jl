using Statistics
using HDF5
using Measurements


"""
    bump_to_nearest_Δx(val, Δx)

Given a value `val`, return `val` rounded to the nearest `Δx`
"""
function bump_to_nearest_Δx(val, Δx)
    return round((val ÷ Δx)*Δx, digits=2)
end


"""
    get_new_bounds(xmin, xmax, ymin, ymax; Δx=0.1)

Given bounding box and a resolution `Δx`, expand bounding box to nearest Δx.
"""
function get_new_bounds(xmin, xmax, ymin, ymax; Δx=0.1)
    # Δx_cm = 100*Δx
    # @assert 100%Δx_cm == 0

    # pad to nearest Δx
    xmin -= Δx
    xmax += Δx
    ymin -= Δx
    ymax += Δx

    return bump_to_nearest_Δx(xmin, Δx), bump_to_nearest_Δx(xmax, Δx), bump_to_nearest_Δx(ymin, Δx), bump_to_nearest_Δx(ymax, Δx)
end



function get_resampled_grid(hsi::HyperspectralImage; Δx=0.1)
    # reinterpolating the results to a new grid
    xmin, xmax = extrema(hsi.X[1,:,:])
    ymin, ymax = extrema(hsi.X[2,:,:])

    bump_to_nearest_Δx(xmin, Δx)

    xmin, xmax, ymin, ymax = get_new_bounds(xmin, xmax, ymin, ymax; Δx=Δx)

    # generate a new x-y grid at the desired resolution
    xs_new = round.(range(xmin, stop=xmax, step=Δx), digits=2)
    ys_new = round.(range(ymin, stop=ymax, step=Δx), digits=2)


    # Xout = permutedims(cat([x for x ∈ xs_new, y ∈ ys_new], [y for x ∈ xs_new, y ∈ ys_new], dims=3), (3,1,2))  # 1.7 ms

    Xhsi = bump_to_nearest_Δx.(hsi.X[1:2,:,:],Δx)

    # generate hsi pixel indices in outbound grid
    Xhsi_is = Matrix{Int}(undef, size(Xhsi, 2), size(Xhsi,3));
    Xhsi_js = Matrix{Int}(undef, size(Xhsi, 2), size(Xhsi,3));

    @tturbo for j ∈ axes(Xhsi, 3), i ∈ axes(Xhsi,2)
        Xhsi_is[i,j] = Int(round((Xhsi[1,i,j] - xmin)/Δx + 1))
        Xhsi_js[i,j] = Int(round((Xhsi[2,i,j] - ymin)/Δx + 1))
    end

    # generate boundary mask
    IsInbounds = [false for i ∈ 1:length(xs_new), j ∈ 1:length(ys_new)];

    @tturbo for j ∈ axes(Xhsi,3), i ∈ axes(Xhsi,2)
        IsInbounds[Xhsi_is[i,j], Xhsi_js[i,j]] = true
    end

    # generate idx mappings
    ij_pixels = findall(IsInbounds)
    ij_notpixels = findall(.!IsInbounds)

    idx_dict = Dict()
    for ij ∈ ij_pixels
        idx_dict[ij] = CartesianIndex{2}[]
    end

    for j ∈ axes(Xhsi_is, 2), i ∈ axes(Xhsi_is,1)
        ij = CartesianIndex(Xhsi_is[i,j], Xhsi_js[i,j])
        push!(idx_dict[ij], CartesianIndex{2}(i,j))
    end


    return xs_new, ys_new, Xhsi_is, Xhsi_js, IsInbounds, ij_pixels, ij_notpixels, idx_dict
end




function resample_datacube_fast(hsi::HyperspectralImage; Δx=0.10)
    # 1. resample to a square grid
    println("\tgenerating new grid")
    xs_new, ys_new, Xhsi_is, Xhsi_js, IsInbounds, ij_pixels, ij_notpixels, idx_dict = get_resampled_grid(hsi;Δx=Δx)

    # 2. allocate latitudes and longitudes
    Longitudes_out= Matrix{Float64}(undef, size(IsInbounds)...)
    Latitudes_out = Matrix{Float64}(undef, size(IsInbounds)...)


    lla_from_utm = LLAfromUTM(hsi.zone, hsi.isnorth, wgs84)

    Threads.@threads for j ∈ axes(Longitudes_out,2)
        for i ∈ axes(Longitudes_out,1)
            lla = lla_from_utm(UTM(xs_new[i], ys_new[j], 0.0))

            Longitudes_out[i, j] = lla.lon
            Latitudes_out[i, j] = lla.lat
        end
    end


    # 3. create vector of labels

    varnames = [
        ["R_"*lpad(idx, 3,"0") for idx ∈ 1:length(hsi.λs)]...,
        "roll",
        "pitch",
        "heading",
        "view_angle",
        "solar_azimuth",
        "solar_elevation",
        "solar_zenith",
    ]

    printnames = [
        ["Reflectance Band "*lpad(idx, 3,"0") for idx ∈ 1:length(hsi.λs)]...,
        "Roll",
        "Pitch",
        "Heading",
        "Viewing Angle",
        "Solar Azimuth",
        "Solar Elevation",
        "Solar Zenith",
    ]

    # 4. Allocate Data Matrices
    Data= Array{Float64}(undef, length(varnames), size(IsInbounds)...);

    # set all out-of-bounds pixels to NaN
    Data[:,ij_notpixels] .= NaN;

    # 5. set up band indices
    ks_reflectance = 1:length(hsi.λs)
    k_roll = findfirst(varnames .== "roll")
    k_pitch = findfirst(varnames .== "pitch")
    k_heading = findfirst(varnames .== "heading")
    k_view = findfirst(varnames .== "view_angle")
    k_az = findfirst(varnames .== "solar_azimuth")
    k_el = findfirst(varnames .== "solar_elevation")
    k_zen = findfirst(varnames .== "solar_zenith")

    # convert Datacube to float first
    Datacube = Float64.(hsi.Datacube)
    # 6. Resample the data
    println("\tinterpolating...")
    Threads.@threads for k ∈ 1:length(ij_pixels)
        @inbounds ij = ij_pixels[k]

        # copy reflectance
        # @inbounds Data[ks_reflectance, ij] = mean(hsi.Datacube[:, idx_dict[ij]], dims=2)
        @inbounds Data[ks_reflectance, ij] = mean(Datacube[:, idx_dict[ij]], dims=2)

        # copy Roll
        @inbounds Data[k_roll, ij] = mean(hsi.Roll[idx_dict[ij]])

        # copy Pitch
        @inbounds Data[k_pitch, ij] = mean(hsi.Pitch[idx_dict[ij]])

        # copy Heading
        @inbounds Data[k_heading, ij] = mean(hsi.Heading[idx_dict[ij]])

        # copy Viewing Angle
        @inbounds Data[k_view, ij] = mean(hsi.ViewAngle[idx_dict[ij]])

        # copy Solar Azimuth
        @inbounds Data[k_az, ij] = mean(hsi.SolarAzimuth[idx_dict[ij]])

        # copy Solar Elevation
        @inbounds Data[k_el, ij] = mean(hsi.SolarElevation[idx_dict[ij]])

        # copy Solar Zenith
        @inbounds Data[k_zen, ij] = mean(hsi.SolarZenith[idx_dict[ij]])
    end


    return xs_new, ys_new, hsi.isnorth, hsi.zone, Longitudes_out, Latitudes_out, IsInbounds, varnames, printnames, hsi.λs, Data

end





function resample_datacube(hsi::HyperspectralImage; Δx=0.10)
    # 1. resample to a square grid
    println("\tgenerating new grid")
    xs_new, ys_new, Xhsi_is, Xhsi_js, IsInbounds, ij_pixels, ij_notpixels, idx_dict = get_resampled_grid(hsi;Δx=Δx)

    # 2. allocate latitudes and longitudes
    Longitudes_out= Matrix{Float64}(undef, size(IsInbounds)...)
    Latitudes_out = Matrix{Float64}(undef, size(IsInbounds)...)


    lla_from_utm = LLAfromUTM(hsi.zone, hsi.isnorth, wgs84)

    Threads.@threads for j ∈ axes(Longitudes_out,2)
        for i ∈ axes(Longitudes_out,1)
            lla = lla_from_utm(UTM(xs_new[i], ys_new[j], 0.0))

            Longitudes_out[i, j] = lla.lon
            Latitudes_out[i, j] = lla.lat
        end
    end


    # 3. create vector of labels
    # 3. create vector of labels

    varnames = [
        ["R_"*lpad(idx, 3,"0") for idx ∈ 1:length(hsi.λs)]...,
        "roll",
        "pitch",
        "heading",
        "view_angle",
        "solar_azimuth",
        "solar_elevation",
        "solar_zenith",
        "mNDWI",
        "NDVI",
        "SR",
        "EVI",
        "AVRI",
        "NDVI_705",
        "MSR_705",
        "MNDVI",
        "VOG1",
        "VOG2",
        "VOG3",
        "PRI",
        "SIPI",
        "PSRI",
        "CRI1",
        "CRI2",
        "ARI1",
        "ARI2",
        "WBI",
        "MCRI",
        "TCARI",
    ]

    printnames = [
        ["Reflectance Band "*lpad(idx, 3,"0") for idx ∈ 1:length(hsi.λs)]...,
        "Roll",
        "Pitch",
        "Heading",
        "Viewing Angle",
        "Solar Azimuth",
        "Solar Elevation",
        "Solar Zenith",
        "Modified Normalized Difference Water Index",
        "Normalized Difference Vegetative Index",
        "Simple Ratio",
        "Enhanced Vegetative Index",
        "Atmospherically Resistant Vegetative Index",
        "Red Edge Normalized Difference Vegetation Index",
        "Modified Red Edge Simple Ratio",
        "Modified Red Edge Normalized Vegetation Index",
        "Vogelmann Red Edge Index",
        "Vogelmann Red Edge Index 2",
        "Vogelmann Red Edge Index 3",
        "Photochemical Reflectance Index",
        "Structure Intensive Pigment Index",
        "Plant Senescence Reflectance Index",
        "Carotenoid Reflectance Index",
        "Carotenoid Reflectance Index 2",
        "Anthocyanin Reflectance Index",
        "Anthocyanin Reflectance Index 2",
        "Water Band Index",
        "Modified Chlorophyll Absorption Reflectance Index",
        "Transformed Chlorophyll Absorption Reflectance Index",
    ]

    # 4. Allocate Data Matrices
    Data= Array{Float64}(undef, length(varnames), size(IsInbounds)...);

    # set all out-of-bounds pixels to NaN
    Data[:,ij_notpixels] .= NaN;

    # 5. set up band indices
    ks_reflectance = 1:length(hsi.λs)
    k_roll = findfirst(varnames .== "roll")
    k_pitch = findfirst(varnames .== "pitch")
    k_heading = findfirst(varnames .== "heading")
    k_view = findfirst(varnames .== "view_angle")
    k_az = findfirst(varnames .== "solar_azimuth")
    k_el = findfirst(varnames .== "solar_elevation")
    k_zen = findfirst(varnames .== "solar_zenith")


    # 6. Resample the data
    println("\tinterpolating...")
    Threads.@threads for k ∈ 1:length(ij_pixels)
        @inbounds ij = ij_pixels[k]

        # copy reflectance
        @inbounds Data[ks_reflectance, ij] = mean(hsi.Datacube[:, idx_dict[ij]], dims=2)[:, 1]

        # copy Roll
        @inbounds Data[k_roll, ij] = mean(hsi.Roll[idx_dict[ij]])

        # copy Pitch
        @inbounds Data[k_pitch, ij] = mean(hsi.Pitch[idx_dict[ij]])

        # copy Heading
        @inbounds Data[k_heading, ij] = mean(hsi.Heading[idx_dict[ij]])

        # copy Viewing Angle
        @inbounds Data[k_view, ij] = mean(hsi.ViewAngle[idx_dict[ij]])

        # copy Solar Azimuth
        @inbounds Data[k_az, ij] = mean(hsi.SolarAzimuth[idx_dict[ij]])

        # copy Solar Elevation
        @inbounds Data[k_el, ij] = mean(hsi.SolarElevation[idx_dict[ij]])

        # copy Solar Zenith
        @inbounds Data[k_zen, ij] = mean(hsi.SolarZenith[idx_dict[ij]])
    end


    return xs_new, ys_new, hsi.isnorth, hsi.zone, Longitudes_out, Latitudes_out, IsInbounds, varnames, printnames, hsi.λs, Data

end




function generate_derived_metrics!(Data, IsInbounds, varnames, λs)
    ij_pixels = findall(IsInbounds)
    k_zen = findfirst(varnames .== "solar_zenith")

    k_λ = k_zen+1

    # mNDWI
    kgreen = findfirst(λs .> 550)
    kswir = findfirst(λs .> 770)

    green = @view Data[kgreen,:,:]
    swir = @view Data[kswir,:,:]

    numer = green .- swir
    denom = green .+ swir

    Data[k_λ,:,:] = numer ./ denom

    # NDVI "normalized difference vegetative index ∈ [-1, 1]"
    k_λ += 1
    kir = findfirst(λs .> 800)
    kred = findfirst(λs .> 680)

    ir_band = @view Data[kir, :,:]
    red_band = @view Data[kred, :,:]

    numer = (ir_band .- red_band)
    denom = (ir_band .+ red_band)

    Data[k_λ,:,:] = numer ./ denom

    # SR "simple ratio ∈ [0, 30]"
    k_λ += 1
    ir_band = @view Data[kir, :,:]
    red_band = @view Data[kred, :,:]

    numer = ir_band
    denom = red_band

    Data[k_λ,:,:] = numer ./ denom

    # EVI "enhanced vegetative index ∈ [-1, 1]"
    k_λ += 1
    kblue = findfirst(λs .> 450)

    ir_band = @view Data[kir, :,:]
    red_band = @view Data[kred, :,:]
    blue_band = @view Data[kblue, :,:]

    numer = 2.5 * (ir_band .- red_band)
    denom = ir_band .+ 6 .* red_band .- 7.5 .* blue_band

    Data[k_λ,:,:] = numer ./ denom

    # AVRI "Atmospherical Resistant Vegitative Indes"
    k_λ += 1

    ir_band = @view Data[kir, :,:]
    red_band = @view Data[kred, :,:]
    blue_band = @view Data[kblue, :,:]

    numer = (ir_band .- 2 .* red_band .+ blue_band)
    denom = (ir_band .+ 2 .* red_band .- blue_band)

    Data[k_λ,:,:] = (numer ./ denom)


    # # NDVI_705 "Red Edge Normalized Difference Vegetation Index"
    k_λ += 1
    kir = findfirst(λs .> 750)
    kred = findfirst(λs .> 705)

    ir_band = @view Data[kir, :,:]
    red_band = @view Data[kred, :,:]

    numer = (ir_band .- red_band)
    denom = (ir_band .+ red_band)

    Data[k_λ,:,:] = (numer ./ denom)


    # # MSR_705 "Modified Red Edge Simple Ratio Index"
    k_λ += 1
    kblue = findfirst(λs .> 445)

    ir_band = @view Data[kir, :,:]
    red_band = @view Data[kred, :,:]
    blue_band = @view Data[kblue, :,:]

    numer = (ir_band .- blue_band)
    denom = (red_band .- blue_band)

    Data[k_λ,:,:] = (numer ./ denom)


    # # MNDVI "modified red edge normalized vegetation index"
    k_λ += 1

    ir_band = @view Data[kir, :,:]
    red_band = @view Data[kred, :,:]
    blue_band = @view Data[kblue, :,:]

    numer = (ir_band .- red_band)
    denom = (ir_band .+ red_band .- 2 .* blue_band)

    Data[k_λ,:,:] =  (numer ./ denom)


    # # VOG1 "vogelmann red edge index"
    k_λ += 1
    kir = findfirst(λs .> 740)
    kred = findfirst(λs .> 720)

    ir_band = @view Data[kir, :,:]
    red_band = @view Data[kred, :,:]

    numer = ir_band
    denom = red_band

    Data[k_λ,:,:] = (numer ./ denom)


    # # VOG2 "vogelmann red edge index 2"
    k_λ += 1
    k1 = findfirst(λs .> 734)
    k2 = findfirst(λs .> 747)
    k3 = findfirst(λs .> 715)
    k4 = findfirst(λs .> 726)

    band1 = @view Data[k1, :,:]
    band2 = @view Data[k2, :,:]
    band3 = @view Data[k3, :,:]
    band4 = @view Data[k4, :,:]

    numer = (band1 .- band2)
    denom = (band3 .+ band4)

    Data[k_λ,:,:] =  (numer ./ denom)


    # # VOG3 "vogelmann red edge index 3 ∈ [0, 20]"
    k_λ += 1
    k1 = findfirst(λs .> 734)
    k2 = findfirst(λs .> 747)
    k3 = findfirst(λs .> 715)
    k4 = findfirst(λs .> 720)

    band1 = @view Data[k1, :,:]
    band2 = @view Data[k2, :,:]
    band3 = @view Data[k3, :,:]
    band4 = @view Data[k4, :,:]

    Data[k_λ,:,:] =  (numer ./ denom)

    # # PRI "photochemical reflectance index" ∈[-1, 1]
    k_λ
    k1 = findfirst(λs .> 531)
    k2 = findfirst(λs .> 570)

    band1 = @view Data[k1, :,:]
    band2 = @view Data[k2, :,:]

    numer = (band1 .- band2)
    denom = (band1 .+ band2)

    Data[k_λ, :,:] = (numer ./ denom)

    # # SIPI "structure intensive pigment index" ∈[0, 2]
    k_λ += 1
    k1 = findfirst(λs .> 800)
    k2 = findfirst(λs .> 445)
    k3 = findfirst(λs .> 680)

    band1 = @view Data[k1, :,:]
    band2 = @view Data[k2, :,:]
    band3 = @view Data[k3, :,:]

    numer = (band1 .- band2)
    denom = (band1 .+ band3)

    Data[k_λ,:,:] = (numer ./ denom)

    # # PSRI "Plant Senescence Reflectance Index"
    k_λ += 1
    k1 = findfirst(λs .> 680)
    k2 = findfirst(λs .> 500)
    k3 = findfirst(λs .> 750)

    band1 = @view Data[k1, :,:]
    band2 = @view Data[k2, :,:]
    band3 = @view Data[k3, :,:]

    numer = (band1 .- band2)
    denom = band3

    Data[k_λ,:,:] = (numer ./ denom)

    # # CRI1 "carotenoid reflectance index"
    k_λ += 1
    k1 = findfirst(λs .> 510)
    k2 = findfirst(λs .> 550)

    band1 = @view Data[k1, :,:]
    band2 = @view Data[k2, :,:]

    Data[k_λ,:,:] = ((1 ./ band1) .- (1 ./ band2))

    # CRI2 "carotenoid reflectance index 2"
    k_λ += 1
    k1 = findfirst(λs .> 510)
    k2 = findfirst(λs .> 700)

    band1 = @view Data[k1, :,:]
    band2 = @view Data[k2, :,:]

    Data[k_λ, :,:] = ((1 ./ band1) .- (1 ./ band2))

    # # ARI1 "anthocyanin reflectance index"
    k_λ += 1
    k1 = findfirst(λs .> 550)
    k2 = findfirst(λs .> 700)

    band1 = @view Data[k1, :,:]
    band2 = @view Data[k2, :,:]

    Data[k_λ,:,:]  = ((1 ./ band1) .- (1 ./ band2))

    # # ARI2 "anthocyanin reflectance index 2"
    k_λ += 1
    k1 = findfirst(λs .> 550)
    k2 = findfirst(λs .> 700)
    k3 = findfirst(λs .> 800)

    band1 = @view Data[k1, :,:]
    band2 = @view Data[k2, :,:]
    band3 = @view Data[k3, :,:]

    Data[k_λ,:,:] = (band3 .* ((1 ./ band1) .- (1 ./ band2)))

    # # WBI "water band index"
    k_λ += 1
    k1 = findfirst(λs .> 900)
    k2 = findfirst(λs .> 970)

    band1 = @view Data[k1, :,:]
    band2 = @view Data[k2, :,:]

    numer = band1
    denom = band2

    Data[k_λ, :,:] = (numer ./ denom)

    # # MCRI "Modified Chlorophyll Absorption Reflectance Index"
    k_λ += 1
    k1 = findfirst(λs .> 550)
    k2 = findfirst(λs .> 670)
    k3 = findfirst(λs .> 701)
    k4 = findfirst(λs .> 780)

    band1 = @view Data[k1, :,:]
    band2 = @view Data[k2, :,:]
    band3 = @view Data[k3, :,:]
    band4 = @view Data[k4, :,:]

    Data[k_λ,:,:] = (((band3 .- band2) .- 0.2 .* (band3 .- band1)) .* (band3 ./ band2))


    # # TCARI "transformed chlorophyll absorption reflectance index"
    k_λ += 1

    band1 = @view Data[k1, :,:]
    band2 = @view Data[k2, :,:]
    band3 = @view Data[k3, :,:]

    Data[k_λ,:,:] = (3 .* ((band3 .- band2) .- 0.2 .* (band3 .- band1) .* (band3 ./ band2)))



end




function save_resampled_hsi(
    xs,
    ys,
    isnorth,
    zone,
    Longitudes,
    Latitudes,
    IsInbounds,
    varnames,
    printnames,
    λs,
    Data,
    outpath;
    Δx=0.10,
    is_spec_chunked=false,
    is_band_chunked=false
    )

    h5open(outpath, "w") do f
        d = create_group(f, "data-Δx_$(Δx)")

        # save resolution
        d["Δx"] = Δx

        # write the UTM coords
        d["X"] = xs
        d["Y"] = ys
        d["isnorth"] = isnorth
        d["zone"] = zone

        # write Lat/Lon grid
        d["Longitudes"] = Longitudes
        d["Latitudes"] = Latitudes

        # write pixel-mask
        d["IsInbounds"] = IsInbounds

        # write variable names
        d["varnames"] = varnames
        d["printnames"] = printnames

        # save wavelength values
        d["λs"] = λs


        # save mean data
        if is_spec_chunked
            d["Data", chunk=(length(varnames),1,1)] = Data
        elseif is_band_chunked
            d["Data", chunk=(1,size(Longitudes)...)] = Data
        else
            d["Data"] = Data
        end


    end
end
