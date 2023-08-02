function Rotation(ϕ,θ,ψ)
    return @SMatrix [
        cos(ψ)*cos(θ)  cos(ψ)*sin(θ)*sin(ϕ)-sin(ψ)*cos(ϕ)  cos(ψ)*sin(θ)*cos(ϕ)+sin(ψ)*sin(ϕ)
        sin(ψ)*cos(θ)  sin(ψ)*sin(θ)*sin(ϕ)+cos(ψ)*cos(ϕ)  sin(ψ)*sin(θ)*cos(ϕ)-cos(ψ)*sin(ϕ)
        -sin(θ)        cos(θ)*sin(ϕ)                       cos(θ)*cos(ϕ)
    ]
end



"""
    generateCoords!(h5::HDF5.File; θ_view = 30.8, z_ground = 292.0, isflipped=false)

Given HSI data stored in an HDF5 file `h5`, georectify the image to produce coordinates for each pixel using position and orientation data from the IMU. This method assumes a flat ground at height `z_ground`. Additionally, the user should specify the viewing angle for their imager in degrees with `θ_view`. Finally, twiddle the witch `isflipped` to deal with potential mix-match
"""
function generateCoords!(h5::HDF5.File; θ_view = 30.8, z_ground = 292.0, isflipped=false)
    nbands = read(h5["raw/radiance/nbands"])  # number of wavelength bins
    samples = read(h5["raw/radiance/ncols"])  # number of pixels per line
    lines = read(h5["raw/radiance/nrows"])    # number of scan lines

    # pre-allocate output arrays
    pixelCoords = Array{Float64}(undef, 3, lines, samples)   # x, y, z in UTMZ
    pixelTimes = Matrix{Float64}(undef, lines, samples)      # time a pixel was captured
    viewingGeom = Array{Float64}(undef, 3, lines, samples)   # store roll, pitch, heading

    pixelLongitudes = Matrix{Float64}(undef, lines, samples) # output longitudes
    pixelLatitudes = Matrix{Float64}(undef, lines, samples)  # output latitudes

    N = samples
    focal_length = ((N-1)/2)/tand(θ_view/2) # compute the focal length in units of "pixels". tand is tangent in degrees

    # loop through and set pixelCoordinates to coordinates in sensor frame
    if !isflipped
        @turbo for k ∈ axes(pixelCoords,1), i ∈ axes(pixelCoords,2), j ∈ axes(pixelCoords,3)
            pixelCoords[k,i,j] = ifelse(k==1, 0, ifelse(k==2, -(N-1)/2 + (j - 1), focal_length))
        end
    else
        @turbo for k ∈ axes(pixelCoords,1), i ∈ axes(pixelCoords,2), j ∈ axes(pixelCoords,3)
            pixelCoords[k,i,j] = ifelse(k==1, 0, ifelse(k==2, (N-1)/2 - (j - 1), focal_length))
        end
    end


    #  read in drone coordinates in  UTM (meters)
    droneCoords = Array{Float64}(undef, 3, lines)
    droneCoords[1,:] .= read(h5["raw/lcf/x"])
    droneCoords[2,:] .= read(h5["raw/lcf/y"])
    droneCoords[3,:] .= read(h5["raw/lcf/z"])

    # read in viewing geometry
    ϕ = read(h5["raw/lcf/roll"])
    θ = read(h5["raw/lcf/pitch"])
    ψ = read(h5["raw/lcf/heading"])

    # read in times
    ts = read(h5["raw/lcf/times"])

    # read other UTMZ information
    isnorth = read(h5["raw/lcf/isnorth"])[1]
    zone = read(h5["raw/lcf/zone"])[1]


    # now we should be able to go through each scanline and compute the updated coordinates...

    pxTmp = copy(pixelCoords[:,1,:])

    Threads.@threads for line ∈ axes(pixelCoords,2)
        # compute scale factor
        #s = (droneCoords[3,line]-z_ground)/focal_length
        s = (droneCoords[3,line]-z_ground)/(focal_length*cos(θ[line]))

        # apply rotation matrix
        mul!(pxTmp, Rotation(ϕ[line], θ[line], ψ[line]), pixelCoords[:,line,:])
        pixelCoords[:,line,:] .= pxTmp

        # pixelCoords[:,line,:] .= Rotation(ϕ[line], θ[line], ψ[line]) .* pixelCoords[:,line,:]

        # scale result to ground scale and shift by drone position
        pixelCoords[:,line,:] .= s .* pixelCoords[:,line,:] .+ droneCoords[:,line]

        # update viewing geometry
        viewingGeom[1,line,:] .= ϕ[line]
        viewingGeom[2,line,:] .= θ[line]
        viewingGeom[3,line,:] .= ψ[line]

        # update pixel times
        pixelTimes[line,:] .= ts[line]

        for j ∈ axes(pixelCoords,3)
            r = LLAfromUTMZ(wgs84)(UTMZ(
                pixelCoords[1,line,j],
                pixelCoords[2,line,j],
                pixelCoords[3,line,j],
                zone,
                isnorth
            ))

            pixelLongitudes[line,j] = r.lon
            pixelLatitudes[line,j] = r.lat
        end
    end

    # now we can add new data to the h5 file
    g = create_group(h5["raw"], "georectified")

    g["X"] = pixelCoords[1,:,:]
    g["Y"] = pixelCoords[2,:,:]
    g["Z"] = pixelCoords[3,:,:]

    g["roll"] = viewingGeom[1,:,:]
    g["pitch"] = viewingGeom[2,:,:]
    g["heading"] = viewingGeom[3,:,:]

    g["times"] = pixelTimes[:,:]

    # generate lat lon grid
    g["longitude"] = pixelLongitudes
    g["latitude"] = pixelLatitudes
end

