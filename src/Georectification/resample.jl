using Statistics
using HDF5
using Trapz

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
        "DVI",
        "EVI",
        "GEMI",
        "GARI",
        "GCI",
        "GDVI",
        "GLI",
        "GNDVI",
        "GOSAVI",
        "GRVI",
        "GSAVI",
        "IPVI",
        "LAI",
        "LCI",
        "MNLI",
        "MSAVI2",
        "Modified Simple Ratio",
        "NLI",
        "NDVI",
        "NPCI",
        "OSAVI",
        "RDVI",
        "SAVI",
        "SR",
        "TDVI",
        "TGI",
        "VARI",
        "WDRVI",
        "AVRI",
        "MCARI",
        "MCARI2",
        "MRENDVI",
        "MRESER",
        "MTVI",
        "RENDVI",
        "TCARI",
        "TVI",
        "VREI1",
        "VREI2",
        "VREI3",
        "PRI",
        "SIPI",
        "SIPI1",
        "PSRI",
        "ARI1",
        "ARI2",
        "CRI1",
        "CRI2",
        "NDWI1",
        "NDWI2",
        "MNDWI",
        "WBI",
        "ACI",
        "MARI",
        "MSI",
        "MTCI",
        "NDII",
        "NDRE",
        "RGRI",
        "RVSI",
        "yaw_minus_azimuth",
        "Σrad",
        "Σdownwelling",
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
        "Difference Vegetation Index",
        "Enhanced Vegetation Index",
        "Global Environmental Monitoring Index",
        "Green Atmospherically Resistant Index",
        "Green Chlorophyll Index",
        "Green Difference Vegetation Index",
        "Green Leaf Index",
        "Green Normalized Difference Vegetation Index",
        "Green Optimized Soil Adjusted Vegetation Index",
        "Green Ratio Vegetation Index",
        "Green Soil Adjusted Vegetation Index",
        "Infrared Percentage Vegetation Index",
        "Leaf Area Index",
        "Leaf Chlorphyll Index",
        "Modified Non-Linear Index",
        "Modified Soil Adjusted Vegetation Index 2",
        "Modified Simple Ratio",
        "Non-Linear Index",
        "Normalized Difference Vegetation Index",
        "Normalized Pigment Chlorophyll Index",
        "Optimized Soil Ajusted Vegetation Index",
        "Renormalized Difference Vegetation Index",
        "Simple Ratio",
        "Transformed Difference Vegetation Index",
        "Triangular Greenness Index",
        "Visible Atmospherically Resistant Index",
        "Wide Dynamic Range Vegetation Index",
        "Atmospherically Resistant Vegetation Index",
        "Modified Chlorophyll Absorption Ratio Index",
        "Modified Chlorophyll Absorption Ratio Index Improved",
        "Modified Red Edge Normalized Difference Vegetation Index",
        "Modified Red Edge Simple Ratio",
        "Modified Triangular Vegetation Index",
        "Red Edge Normalized Difference Vegetation Index",
        "Transformed Chlorophyll Absorption Reflectance Index",
        "Triangular Vegetation Index",
        "Vogelmann Red Edge Index 1",
        "Vogelmann Red Edge Index 2",
        "Vogelmann Red Edge Index 3",
        "Photochemical Reflectance Index",
        "Structure Insensitive Pigment Index",
        "Structure Independent Pigment Index",
        "Plant Senescence Reflectance Index",
        "Anthocyanin Reflectance Index 1",
        "Anthocyanin Reflectance Index 2",
        "Carotenoid Reflectance Index 1",
        "Carotenoid Reflectance Index 2",
        "Normalized Difference Water Index 1",
        "Normalized Difference Water Index 2",
        "Modified Normalized Difference Water Index",
        "Water Band Index",
        "Anthocyanin Content Index",
        "Modified Anthocyanin Reflectance Index",
        "Moisture Stress Index",
        "MERIS Terrestrial Chlorophyll Index",
        "Normalized Difference Infrared Index",
        "Normalized Difference Red Edge",
        "Red Green Ratio Index",
        "Red Edge Vegetation Stress Index",
        "Heading - Solar Azimuth",
        "Total Pixel Intensity",
        "Total Downwelling Intensity",
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
    k_yma = findfirst(varnames .== "yaw_minus_azimuth")
    k_tot = findfirst(varnames .== "Σrad")
    k_down = findfirst(varnames .== "Σdownwelling")

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
        @inbounds Data[k_az, ij] = mean(hsi.SolarAzimuth[idx_dict[ij]]) *  π / 180.0

        # copy Solar Elevation
        @inbounds Data[k_el, ij] = mean(hsi.SolarElevation[idx_dict[ij]])  *  π / 180.0

        # copy Solar Zenith
        @inbounds Data[k_zen, ij] = mean(hsi.SolarZenith[idx_dict[ij]])  *  π / 180.0

        # compute heading difference
        @inbounds Data[k_yma, ij] = Data[k_heading, ij] - Data[k_az, ij]

        # integrate to obtain per pixel radiance in W/m^2
        @inbounds Data[k_tot, ij] = trapz(hsi.λs .* 1e-3,  π .* Data[ks_reflectance, ij] .* 1e-2)

        # set initial value of downwelling to zero
        @inbounds Data[k_down, ij] = 0.0
    end


    return xs_new, ys_new, hsi.isnorth, hsi.zone, Longitudes_out, Latitudes_out, IsInbounds, varnames, printnames, hsi.λs, Data

end




function generate_derived_metrics!(Data, IsInbounds, varnames, λs)
    ij_pixels = findall(IsInbounds)
    k_zen = findfirst(varnames .== "solar_zenith")



    # use ENVI definitions here: https://www.nv5geospatialsoftware.com/docs/spectralindices.html
    # we choose the band "center"
    #
    # band  | minimum | center  | maximum
    # -----------------------------------
    # blue  | 400 nm  | 470 nm  | 500 nm
    # green | 500 nm  | 550 nm  | 600 nm
    # red   | 600 nm  | 650 nm  | 700 nm
    # NIR   | 760 nm  | 860 nm  | 960 nm
    # SWIR1 | 1550 nm | 1650 nm | 1750 nm

    λ_blue = 470
    λ_green = 550
    λ_red = 650


    kblue = findfirst(λs .>= 770)
    kgreen = findfirst(λs .>= 550)
    kred = findfirst(λs .>= 650)
    knir = findfirst(λs .>= 860)
    kswir = findfirst(λs .>= 1009)  # an *almost* shortwave

    k_430 = findfirst(λs .>= 430)
    k_445 = findfirst(λs .>= 445)
    k_450 = findfirst(λs .>= 450)
    k_500 = findfirst(λs .>= 500)
    k_510 = findfirst(λs .>= 510)
    k_531 = findfirst(λs .>= 531)
    k_550 = findfirst(λs .>= 550)
    k_570 = findfirst(λs .>= 570)
    k_670 = findfirst(λs .>= 670)
    k_680 = findfirst(λs .>= 680)
    k_700 = findfirst(λs .>= 700)
    k_705 = findfirst(λs .>= 705)
    k_710 = findfirst(λs .>= 710)
    k_714 = findfirst(λs .>= 714)
    k_715 = findfirst(λs .>= 715)
    k_720 = findfirst(λs .>= 720)
    k_726 = findfirst(λs .>= 726)
    k_733 = findfirst(λs .>= 733)
    k_734 = findfirst(λs .>= 734)
    k_740 = findfirst(λs .>= 740)
    k_747 = findfirst(λs .>= 747)
    k_750 = findfirst(λs .>= 750)
    k_752 = findfirst(λs .>= 752)
    k_790 = findfirst(λs .>= 790)
    k_800 = findfirst(λs .>= 800)
    k_850 = findfirst(λs .>= 850)
    k_900 = findfirst(λs .>= 900)
    k_970 = findfirst(λs .>= 970)


    band_blue = @view Data[kblue,:,:]
    band_green = @view Data[kgreen,:,:]
    band_red = @view Data[kred,:,:]
    band_nir = @view Data[knir,:,:]
    band_swir = @view Data[kswir,:,:]

    band_430 = @view Data[k_430,:,:]
    band_445 = @view Data[k_445,:,:]
    band_450 = @view Data[k_450,:,:]
    band_500 = @view Data[k_500,:,:]
    band_510 = @view Data[k_510,:,:]
    band_531 = @view Data[k_531,:,:]
    band_550 = @view Data[k_550,:,:]
    band_570 = @view Data[k_570,:,:]
    band_670 = @view Data[k_670,:,:]
    band_680 = @view Data[k_680,:,:]
    band_700 = @view Data[k_700,:,:]
    band_705 = @view Data[k_705,:,:]
    band_710 = @view Data[k_710,:,:]
    band_714 = @view Data[k_714,:,:]
    band_715 = @view Data[k_715,:,:]
    band_720 = @view Data[k_720,:,:]
    band_726 = @view Data[k_726,:,:]
    band_733 = @view Data[k_733,:,:]
    band_734 = @view Data[k_734,:,:]
    band_740 = @view Data[k_740,:,:]
    band_747 = @view Data[k_747,:,:]
    band_750 = @view Data[k_750,:,:]
    band_752 = @view Data[k_752,:,:]
    band_790 = @view Data[k_790,:,:]
    band_800 = @view Data[k_800,:,:]
    band_850 = @view Data[k_850,:,:]
    band_900 = @view Data[k_900,:,:]
    band_970 = @view Data[k_970,:,:]

    # list is here: https://www.nv5geospatialsoftware.com/docs/VegetationIndices.html

    # ---------------------------------------
    # Broadband Greenness
    # https://www.nv5geospatialsoftware.com/docs/broadbandgreenness.html
    # ---------------------------------------

    # DVI: Difference Vegetation Index
    k_wav = k_zen+1
    Data[k_wav,:,:] = band_nir .- band_red

    # EVI: Enhanced Vegetation Index
    k_wav += 1

    numer = 2.5 .* (band_nir .- band_red)
    denom = band_nir .+ (6.0 .* band_red) .- (7.5 .* band_blue) .+ 1.0
    Data[k_wav,:,:] = numer ./ denom


    # GEMI: Global Environmental Monitoring Index
    k_wav += 1

    eta = 2 .* (band_nir .^2 .- band_red .^2) .+ (1.5 .* band_nir) .+ (0.5 .* band_red)
    Data[k_wav, :,:] = eta .* (1.0 .- 0.25 .* eta) .- ((band_red .- 1.125) ./ (1.0 .- band_red))


    # GARI: Green Atmospherically Resistant Index
    k_wav += 1
    γ = 1.7
    numer = band_nir .- (band_green .-  γ .* (band_blue .- band_red))
    denom = band_nir .+ (band_green .-  γ .* (band_blue .- band_red))
    Data[k_wav, :,:] = numer ./ denom

    # GCI: green chlorophyll index
    k_wav += 1

    Data[k_wav,:,:] = (band_nir ./ band_green) .- 1.0


    # GDVI: green difference vegetation index
    k_wav += 1

    Data[k_wav,:,:] = band_nir .- band_green


    # GLI: Green Leaf Index
    k_wav += 1
    numer = (band_green .- band_red) .+ (band_green .- band_blue)
    denom = (2.0 .* band_green) .+ band_red .+ band_blue

    Data[k_wav, :,:] = numer ./ denom


    # GNDVI: Green Normalized Difference Vegetation Index
    k_wav += 1
    numer = band_nir .- band_green
    denom = band_nir .+ band_green

    Data[k_wav,:,:] = numer ./ denom

    # GOSAVI: Green Optimized Soil Adjusted Vegetation Index
    k_wav += 1
    numer = band_nir .- band_green
    denom = band_nir .+ band_green .+ 0.16

    Data[k_wav,:,:] = numer./denom

    # GRVI: Green Ratio Vegetation Index
    k_wav += 1
    Data[k_wav,:,:] = band_nir ./ band_green


    # GSAVI: Green soil adjusted vegetation index
    k_wav += 1
    numer = 1.5 .* (band_nir .- band_green)
    denom = band_nir .+ band_green .+ 0.5

    Data[k_wav,:,:] = numer ./ denom


    # IPVI: Infrared Percentage Vegetation Index
    k_wav += 1
    Data[k_wav,:,:] = band_nir ./ (band_nir .+ band_red)


    # LAI: Leaf Area Index
    k_wav += 1
    numer = 2.5 .* (band_nir .- band_red)
    denom = band_nir .+ (6.0 .* band_red) .- (7.5 .* band_blue) .+ 1.0
    Data[k_wav,:,:] = 3.618 .* (numer ./ denom) .- 0.118

    # LCI: Leaf Chlorophyll Index
    k_wav += 1

    numer = band_850 .- band_710
    denom = band_850 .+ band_680
    Data[k_wav,:,:] = numer ./ denom

    # MNLI: Modified Non-Linear Index
    k_wav += 1

    L = 0.5
    numer = (band_nir .^ 2 .- band_red) .* (1 .+ L)
    denom = band_nir .^ 2 .+ band_red .+ L
    Data[k_wav,:,:] = numer ./ denom

    # MSAVI2: Modified Soil Adjusted Vegetation Index 2
    k_wav += 1

    Data[k_wav,:,:] = 0.5 .* (2.0 .* band_nir .+ 1.0 .- sqrt.((2.0 .* band_nir .+ 1.0).^2 .- 8.0 .* (band_nir .- band_red)))


    # MSR: Modified Simple Ratio
    k_wav += 1

    numer = (band_nir ./ band_red) .- 1.0
    denom = sqrt.(band_nir ./ band_red) .+ 1.0

    # NLI: Nonlinear Index
    k_wav += 1

    numer = band_nir .^ 2 .- band_red
    denom = band_nir .^ 2 .+ band_red
    Data[k_wav,:,:] = numer ./ denom

    # NDVI: Normalized Difference Vegetation Index
    k_wav += 1

    numer = band_nir .- band_red
    denom = band_nir .+ band_red
    Data[k_wav,:,:] = numer ./ denom


    # NPCI: Normalized Pigment Chlorophyll Index
    k_wav += 1

    numer = band_680 .- band_430
    denom = band_680 .+ band_430

    Data[k_wav,:,:] = numer ./ denom


    # OSAVI: Optimized Soil Adjusted Vegetation Index
    k_wav += 1

    numer = band_nir .- band_red
    denom = band_nir .+ band_red .+ 0.16

    Data[k_wav,:,:] = numer ./ denom


    # RDVI: renormalized difference vegetation index
    k_wav += 1

    numer = band_nir .- band_red
    denom = sqrt.(band_nir .+ band_red)

    Data[k_wav,:,:] = numer ./ denom


    # SAVI: Soil Adjusted Vegetation Index
    k_wav += 1

    L = 0.5
    numer = 1.5 .* (band_nir .- band_red)
    denom = band_nir .+ band_red .+ 0.5

    Data[k_wav,:,:] = numer ./ denom


    # SR: Simple Ratio
    k_wav += 1

    Data[k_wav,:,:] = band_nir ./ band_red


    # TDVI: Transformed Difference Vegetation Index
    k_wav += 1

    numer = 1.5 .* band_nir .- band_red
    denom = sqrt.(band_nir .^2 .+ band_red .+ 0.5)

    Data[k_wav,:,:] = numer ./ denom


    # TGI: Triangular Greenness Index
    k_wav += 1

    Data[k_wav,:,:] = 0.5 .* ((λ_red- λ_blue) .* (band_red .- band_green) .- (λ_red- λ_green) .* (band_red .- band_blue))


    # VARI: Visible Atmospherically Resistant Index
    k_wav += 1

    numer = band_green .- band_red
    denom = band_green .+ band_red .- band_blue

    Data[k_wav,:,:] = numer ./ denom


    # WDRVI: Wide Dynamic Range Vegetation Index
    k_wav += 1

    a = 0.2
    numer = a .* band_nir .- band_red
    denom = a .* band_nir .+ band_red

    Data[k_wav,:,:] = numer ./ denom


    # ----------------------------
    # Narrow Band Greenness
    # https://www.nv5geospatialsoftware.com/docs/narrowbandgreenness.html
    # ----------------------------

    # ARVI: Atmospherically Resistant Vegetation Index
    k_wav += 1

    γ = 1.0
    numer = band_800 .- (band_680 .-  γ .* (band_450 .- band_680))
    denom = band_800 .+ (band_680 .-  γ .* (band_450 .- band_680))

    Data[k_wav,:,:] = numer ./ denom

    # MCARI: Modified Chlorophyll Absorption Ratio Index
    k_wav += 1

    Data[k_wav,:,:] = ((band_700 .- band_670) .- 2.0 .*(band_700 .- band_550)) .* (band_700 ./ band_670)


    # MCARI2: Modified Chlorophyll Absorption Ratio Index Improved
    k_wav += 1

    numer = 1.5 .* (2.5 .* (band_800 .- band_670) .- 1.3 .* (band_800 .- band_550))
    denom = sqrt.((2.0 .* band_800 .+ 1.0) .^ 2 .- (6.0 .* band_800 .- 5 .* sqrt.(band_670)) .- 0.5)
    Data[k_wav,:,:] = numer ./ denom


    # MRENDVI: Modified Red Edge Normalized Difference Vegetation Index
    k_wav += 1

    numer = band_750 .- band_705
    denom = band_750 .+ band_705 .- (2.0 .* band_445)
    Data[k_wav,:,:] = numer ./ denom


    # MRESER: Modified Red Edge Simple Ratio
    k_wav += 1

    numer = band_750 .- band_445
    denom = band_705 .- band_445
    Data[k_wav,:,:] = numer ./ denom


    # MTVI: Modified Triangular Vegetation Index
    k_wav += 1

    Data[k_wav,:,:] = 1.2 .* (1.2 .* (band_800 .- band_550) - 2.5 .* (band_670 .- band_550))



    # RENDVI: Red Edge Normalized Difference Vegetation Index
    k_wav += 1

    numer = band_750 .- band_705
    denom = band_750 .+ band_705
    Data[k_wav,:,:] = numer ./ denom

    # TCARI: Transformed Chlorophyll Absorption Reflectance Index
    k_wav += 1

    Data[k_wav,:,:] = 3.0 .* ((band_700 .- band_670) .- 0.2 .* (band_700 .- band_550).*(band_700 ./ band_670))


    # TVI: Triangular Vegetation Index
    k_wav += 1

    Data[k_wav,:,:] = 0.5 .* (120 .* (band_750 .- band_550) .- 200 .* (band_670 .- band_550))


    # VREI1: Vogelmann Red Edge Index 1
    k_wav += 1

    Data[k_wav,:,:] = band_740 ./ band_720


    # VREI2: Vogelmann Red Edge Index 2
    k_wav += 1

    numer = (band_734 .- band_747)
    denom = (band_715 .+ band_726)
    Data[k_wav,:,:] =  numer ./ denom


    # VREI3: Vogelmann Red Edge Index 3
    k_wav += 1

    numer = (band_734 .- band_747)
    denom = (band_715 .+ band_720)
    Data[k_wav,:,:] = numer ./ denom


    # ----------------------------
    # Light Use Efficiency
    # https://www.nv5geospatialsoftware.com/docs/lightuseefficiency.html
    # ----------------------------

    # PRI: Photochemical Reflectance Index
    k_wav += 1

    numer = band_531 .- band_570
    denom = band_531 .+ band_570
    Data[k_wav,:,:] = numer ./ denom


    # SIPI: Structure Insensitive Pigment Index
    k_wav += 1

    numer = band_800 .- band_445
    denom = band_800 .+ band_680
    Data[k_wav,:,:] = numer ./ denom

    # SIPI1: Structure Independent Pigment Index
    k_wav += 1

    numer = band_445 .- band_800
    denom = band_670 .- band_800
    Data[k_wav,:,:] = numer ./ denom

    # ----------------------------
    # Dry Senescent Carbon
    # https://www.nv5geospatialsoftware.com/docs/drysenescentcarbon.html
    # ----------------------------

    # PSRI: Plant Senescence Reflectance Index
    k_wav += 1

    numer = band_680 .- band_500
    Data[k_wav,:,:] = numer ./ band_750

    # ----------------------------
    # Leaf Pigments
    # https://www.nv5geospatialsoftware.com/docs/leafpigments.html
    # ----------------------------

    # ARI1: Anthocyanin Reflectance Index 1
    k_wav += 1

    Data[k_wav,:,:] = (1 ./ band_550) .- (1 ./ band_700)


    # ARI2: Anthocyanin Reflectance Index 2
    k_wav += 1

    Data[k_wav,:,:] = band_800 .* ((1 ./ band_550) .- (1 ./ band_700))

    # CRI1: Carotenoid Reflectance Index 1
    k_wav += 1

    Data[k_wav,:,:] = (1 ./ band_510) .- (1 ./ band_550)

    # CRI2: Carotenoid Reflectance Index 2
    k_wav += 1

    Data[k_wav,:,:] = (1 ./ band_510) .- (1 ./ band_700)



    # ----------------------------
    # Canopy Water Content
    # https://www.nv5geospatialsoftware.com/docs/canopywatercontent.html
    # ----------------------------

    # NDWI1: Normalized Difference Water Index 1
    k_wav += 1

    numer = band_green .- band_nir
    denom = band_green .+ band_nir
    Data[k_wav,:,:] = numer ./ denom

    # NDWI2: Normalized Difference Water Index 2
    k_wav += 1

    numer = band_nir .- band_swir
    denom = band_nir .+ band_swir
    Data[k_wav,:,:] = numer ./ denom


    # MNDWI: Modified Normalized Difference Water Index
    k_wav += 1

    numer = band_green .- band_swir
    denom = band_green .+ band_swir
    Data[k_wav,:,:] = numer ./ denom


    # WBI: Water Band Index
    k_wav += 1

    Data[k_wav,:,:] = band_970 ./ band_900




    # ----------------------------
    # Additional Indices from the HSI Book
    # ----------------------------

    # ACI: Anthocyanin Content Index
    k_wav += 1

    Data[k_wav,:,:] = band_green ./ band_nir


    # CIre: Chlorophyll Index Red Edge
    k_wav += 1

    Data[k_wav,:,:] = (band_nir ./ band_705) .- 1.0


    # MARI: Modified Anthocyanin Reflectance Index
    k_wav += 1

    Data[k_wav,:,:] = ((1 ./ band_550) .- (1 ./ band_700)) .* band_nir


    # MSI: Moisture Stress Index
    k_wav += 1

    Data[k_wav,:,:] = band_swir ./ band_nir


    # MTCI: MERIS Terrestrial Chlorophyll Index
    k_wav += 1

    k_1 = findfirst(λs .> 753.75)
    k_2 = findfirst(λs .> 708.75)
    k_3 = findfirst(λs .> 681.25)
    
    band_1 = @view Data[k_1,:,:]
    band_2 = @view Data[k_2,:,:]
    band_3 = @view Data[k_3,:,:]

    Data[k_wav,:,:] = (band_1 .- band_2) ./ (band_2 .- band_3)


    # NDII: Normalized Difference Infrared Index
    k_wav += 1

    Data[k_wav,:,:] = (band_nir .- band_swir) ./ (band_nir .+ band_swir)


    # NDRE: Normalized Difference Red Edge
    k_wav += 1

    numer = band_790 .- band_720
    denom = band_790 .+ band_720
    Data[k_wav,:,:] = numer ./ denom


    # RGRI: Red Green Ratio Index
    k_wav += 1

    Data[k_wav,:,:] = band_red ./ band_green


    # RVSI: Red Edge Vegetation Stress Index
    k_wav += 1

    Data[k_wav,:,:] = (band_714 .+ band_752) ./ 2.0 .- band_733

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
