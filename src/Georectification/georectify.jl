function Rotation(ϕ,θ,ψ)
    return @SMatrix [
        cos(ψ)*cos(θ)  cos(ψ)*sin(θ)*sin(ϕ)-sin(ψ)*cos(ϕ)  cos(ψ)*sin(θ)*cos(ϕ)+sin(ψ)*sin(ϕ)
        sin(ψ)*cos(θ)  sin(ψ)*sin(θ)*sin(ϕ)+cos(ψ)*cos(ϕ)  sin(ψ)*sin(θ)*cos(ϕ)-cos(ψ)*sin(ϕ)
        -sin(θ)        cos(θ)*sin(ϕ)                       cos(θ)*cos(ϕ)
    ]
end


"""
    imagecoords(i, j, N, f)

Given pixel indices `i` and `j`, compute the coordinates of a HSI pixel in the image coordinate system. The height is the focal length `f`, and the pixels are assumed to lie along the y-axis.
"""
function imagecoords(i, j, N, f)
    if i == 1
        return 0.0  # x position
    elseif i == 2
        return (N-1)/2 - (j-1)  # y position
    else
        return f  # z position
    end
end

"""
    imagecoordsFlipped(i, j, N, f)

Same as `imagecoord()` but with the y axis flipped. This is useful in case the camera settings for direction of flight are backwards.
"""
function imagecoordsFlipped(i, j, N, f)
    if i == 1
        return 0.0 # x position
    elseif i == 2
        return -(N-1)/2 + (j-1)  # y position (flipped)
    else
        return f
    end
end



T_n_E = @SMatrix [0 1 0; 1 0 0; 0 0 -1 ]  # matrix to convert navigation system (IMU) to the Earth System (UTM) from Baumker paper page 7



# function generateCoords!(hsi::HyperspectralImage,
#                          fdata::FlightData;
#                          θ_view=30.8,
#                          z_ground=292.0,
#                          isflipped=false
#                          )
#     # update the hsi with the correct start time
#     hsi.start_time[1] = fdata.start_time

#     # we already preallocated the output matrices in the struct
#     bands = nbands(hsi)
#     N = nsamples(hsi)
#     lines = nscanlines(hsi)

#     # compute focal length
#     focal_length = ((N-1)/2)/tand(θ_view/2) # in units of "pixels". tand is tangent in degrees

#     # NOTE: using SMatrix for rs_pixel_sensor isn't great as the matrix is
#     # way too big to be allocated on the stack and still run fast

#     if isflipped
#         # rs_pixel_sensor = SMatrix{3,N}([imagecoordsFlipped(i,j,N,f) for i∈1:3, j∈1:N])
#         rs_pixel_sensor = [imagecoordsFlipped(i,j,N,focal_length) for i∈1:3, j∈1:N]
#     else
#         # rs_pixel_sensor = SMatrix{3,N}([imagecoords(i,j,N,f) for i∈1:3, j∈1:N])
#         rs_pixel_sensor = [imagecoords(i,j,N,focal_length) for i∈1:3, j∈1:N]
#     end


#     # loop through and set viewing angle
#     @tturbo for i ∈ 1:N, j∈1:lines
#         hsi.ViewAngle[i,j] = atan(rs_pixel_sensor[2,i], rs_pixel_sensor[3,i])
#     end

#     # read in viewing geometry
#     ϕ = @view fdata.rolls[:]
#     θ = @view fdata.pitches[:]
#     ψ = @view fdata.headings[:]

#     # read in times
#     ts = @view fdata.times[:]

#     # read other UTMZ information
#     hsi.isnorth[1] = fdata.isnorth
#     hsi.zone[1] = fdata.zone

#     # set up drone coordinates
#     droneCoords = Array{Float64}(undef, 3, lines)
#     droneCoords[1,:] .= fdata.xs
#     droneCoords[2,:] .= fdata.ys
#     droneCoords[3,:] .= fdata.zs

#     Threads.@threads for line ∈ 1:lines
#         # compute scale factor
#         #s = (droneCoords[3,line]-z_ground)/focal_length
#         s = (droneCoords[3,line]-z_ground)/(focal_length*cos(θ[line]))

#         # compute object coords in UTM
#         # T_n_E is conversion from navigation frame to earth frame
#         # Rotation(ϕ,θ,ψ) is conversion from sensor frame to navigation frame
#         # note difference between active/passive transformation
#         rs_object_utm = droneCoords[:,line] .+ s .* T_n_E*Rotation(ϕ[line], θ[line],ψ[line])*rs_pixel_sensor

#         hsi.X[:,:,line] .= rs_object_utm

#         # update viewing geometry
#         hsi.Roll[:,line] .= ϕ[line]
#         hsi.Pitch[:,line] .= θ[line]
#         hsi.Heading[:,line] .= ψ[line]

#         # update pixel times
#         hsi.Times[:,line] .= ts[line]
#     end
# end



function generateCoords!(
    X,
    Longitudes,
    Latitudes,
    Roll,
    Pitch,
    Heading,
    Times,
    ViewAngle,
    fdata::FlightData;
    θ_view=30.8,
    z_ground=292.0,
    isflipped=false
    )


    # we already preallocated the output matrices in the struct
    N = size(X,2)
    lines = size(X,3)

    # compute focal length
    focal_length = ((N-1)/2)/tand(θ_view/2) # in units of "pixels". tand is tangent in degrees

    # NOTE: using SMatrix for rs_pixel_sensor isn't great as the matrix is
    # way too big to be allocated on the stack and still run fast
    if isflipped
        # rs_pixel_sensor = SMatrix{3,N}([imagecoordsFlipped(i,j,N,f) for i∈1:3, j∈1:N])
        rs_pixel_sensor = [imagecoordsFlipped(i,j,N,focal_length) for i∈1:3, j∈1:N]
    else
        # rs_pixel_sensor = SMatrix{3,N}([imagecoords(i,j,N,f) for i∈1:3, j∈1:N])
        rs_pixel_sensor = [imagecoords(i,j,N,focal_length) for i∈1:3, j∈1:N]
    end


    # loop through and set viewing angle
    @tturbo for i ∈ 1:N, j∈1:lines
        ViewAngle[i,j] = atan(rs_pixel_sensor[2,i], rs_pixel_sensor[3,i])
    end

    # read in viewing geometry
    ϕ = @view fdata.rolls[:]
    θ = @view fdata.pitches[:]
    ψ = @view fdata.headings[:]

    # read in times
    ts = @view fdata.times[:]

    # read other UTMZ information
    isnorth = fdata.isnorth
    zone = fdata.zone

    # set up drone coordinates
    droneCoords = Array{Float64}(undef, 3, lines)
    @inbounds droneCoords[1,:] .= fdata.xs
    @inbounds droneCoords[2,:] .= fdata.ys
    @inbounds droneCoords[3,:] .= fdata.zs

    Threads.@threads for line ∈ 1:lines
        # compute scale factor
        #s = (droneCoords[3,line]-z_ground)/focal_length
        s = (droneCoords[3,line]-z_ground)/(focal_length*cos(θ[line]))

        # compute object coords in UTM
        # T_n_E is conversion from navigation frame to earth frame
        # Rotation(ϕ,θ,ψ) is conversion from sensor frame to navigation frame
        # note difference between active/passive transformation
        rs_object_utm = droneCoords[:,line] .+ s .* T_n_E*Rotation(ϕ[line], θ[line],ψ[line])*rs_pixel_sensor

        X[:,:,line] .= rs_object_utm

        # convert back to latitude, longitude
        rs_object_lla = [LLAfromUTMZ(wgs84)(UTMZ(rs_object_utm[:,j]..., zone, isnorth)) for j∈1:N]
        @inbounds Latitudes[:,line] .= [lla.lat for lla ∈ rs_object_lla]
        @inbounds Longitudes[:,line] .= [lla.lon for lla ∈ rs_object_lla]




        # update viewing geometry
        Roll[:,line] .= ϕ[line]
        Pitch[:,line] .= θ[line]
        Heading[:,line] .= ψ[line]

        # update pixel times
        Times[:,line] .= ts[line]
    end
end



