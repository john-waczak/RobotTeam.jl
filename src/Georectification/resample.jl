

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
    Δx_cm = 100*Δx
    @assert 100%Δx_cm == 0

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


    xmin, xmax, ymin, ymax = get_new_bounds(xmin, xmax, ymin, ymax)


    # estimate the bound for minimum pixel spacing in meters
    npix = max(size(hsi.X)...)
    Δx_min = min((ymax-ymin)/npix, (xmax-xmin)/npix)


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

    # generate array w/ number of HSI pixels per location
    Npixels = zeros(Int, size(IsInbounds)...)

    Threads.@threads for j ∈ axes(Xhsi_is,2)
        for i ∈ axes(Xhsi_is, 1)
            Npixels[Xhsi_is[i,j], Xhsi_js[i,j]] += 1
        end
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


    return xs_new, ys_new, Xhsi_is, Xhsi_js, IsInbounds, Npixels, ij_pixels, ij_notpixels, idx_dict
end


# # create vector of labels
# varnames = []
# printnames = []
# for i ∈ 1:length(hsi.λs)
#     idx = lpad(i, 3, "0")
#     push!(varnames, "R_$(idx)")
#     push!(printnames, "Reflectance Band $(idx)")
# end

# push!(varnames, "roll")
# push!(printnames, "Roll")

# push!(varnames, "pitch")
# push!(printnames, "Pitch")

# push!(varnames, "heading")
# push!(printnames, "Heading")

# push!(varnames, "camera_angle")
# push!(printnames, "Camera Angle")

# push!(varnames, "solar_azimuth")
# push!(printnames, "Solar Azimuth")

# push!(varnames, "solar_elevation")
# push!(printnames, "Solar Elevation")

# push!(varnames, "solar_zenith")
# push!(printnames, "Solar Zenith")

# push!(varnames, "mNDWI")
# push!(printnames, "mNDWI")

# push!(varnames, "NDVI")
# push!(printnames, "NDVI")

# push!(varnames, "SR")
# push!(printnames, "SR")

# push!(varnames, "EVI")
# push!(printnames, "EVI")

# push!(varnames, "AVRI")
# push!(printnames, "AVRI")

# push!(varnames, "NDVI_705")
# push!(printnames, "NDVI_705")

# push!(varnames, "MSR_705")
# push!(printnames, "MSR_705")

# push!(varnames, "MNDVI")
# push!(printnames, "MNDVI")

# push!(varnames, "VOG1")
# push!(printnames, "VOG1")

# push!(varnames, "VOG2")
# push!(printnames, "VOG2")

# push!(varnames, "VOG3")
# push!(printnames, "VOG3")

# push!(varnames, "PRI")
# push!(printnames, "PRI")

# push!(varnames, "SIPI")
# push!(printnames, "SIPI")

# push!(varnames, "PSRI")
# push!(printnames, "PSRI")

# push!(varnames, "CRI1")
# push!(printnames, "CRI1")

# push!(varnames, "CRI2")
# push!(printnames, "CRI2")

# push!(varnames, "ARI1")
# push!(printnames, "ARI1")

# push!(varnames, "ARI2")
# push!(printnames, "ARI2")

# push!(varnames, "WBI")
# push!(printnames, "WBI")

# push!(varnames, "MCRI")
# push!(printnames, "MCRI")

# push!(varnames, "TCARI")
# push!(printnames, "TCARI")

