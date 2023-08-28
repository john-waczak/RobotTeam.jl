module Visualization

using ..ENVI: HyperspectralImage

using Images
using HDF5
using GLMakie
using GeometryBasics

export process_image
export generateRGB!
export getRGB
export vis_cube


function process_image(img;α=10.0,β=0.0)
    # see example: https://www.satmapper.hu/en/rgb-images/
    img_out = copy(img)
    # 1. normalize each band
    for i ∈ 1:3
        #brighten
        img_out[i,:,:] .= clamp.(α.*img_out[i,:,:] .+ β, 0, 1)

        # normalize
        # bmin, bmax = extrema(img_out[i,:,:])
        # img_out[i,:,:] .= (img_out[i,:,:] .- bmin)./(bmax - bmin)

    end
    return img_out
end


function getRGB(hsi::HyperspectralImage; λred=630.0, λgreen=532.0, λblue=465.0)
    λs = hsi.λs
    idx_r = argmin(abs.(λs .- λred))
    idx_g = argmin(abs.(λs .- λgreen))
    idx_b = argmin(abs.(λs .- λblue))

    img = hsi.Reflectance[[idx_r, idx_g, idx_b], :, :]

    imgp = process_image(img)

    return colorview(RGB, imgp)
end

function getRGB(Data, λs, ij_pixels; λred=630.0, λgreen=532.0, λblue=465.0)
    img = zeros(4, size(Data,2), size(Data,3))

    idx_r = argmin(abs.(λs .- λred))
    idx_g = argmin(abs.(λs .- λgreen))
    idx_b = argmin(abs.(λs .- λblue))

    Threads.@threads for ij ∈ ij_pixels
        img[1, ij] = Data[idx_r, ij]
        img[2, ij] = Data[idx_g, ij]
        img[3, ij] = Data[idx_b, ij]
        img[4, ij] = 1.0
    end

    imgp = process_image(img)

    return colorview(RGBA, imgp)
end


function getRGB(h5::HDF5.File; λred=630.0, λgreen=532.0, λblue=465.0, Δx=0.10, α=10.0, β=0.0)
    λs = h5["data-Δx_$(Δx)/λs"][:]

    λred=630.0
    λgreen=532.0
    λblue=465.0

    idx_r = argmin(abs.(λs .- λred))
    idx_g = argmin(abs.(λs .- λgreen))
    idx_b = argmin(abs.(λs .- λblue))

    Rr = h5["data-Δx_$(Δx)/Data_μ"][idx_r, :, :]
    Rg = h5["data-Δx_$(Δx)/Data_μ"][idx_g, :, :]
    Rb = h5["data-Δx_$(Δx)/Data_μ"][idx_b, :, :]

    ij_pixels = findall(h5["data-Δx_$(Δx)/IsInbounds"][:,:])
    img = zeros(4, size(Rr)...)

    Threads.@threads for ij ∈ ij_pixels
        img[1, ij] = Rr[ij]
        img[2, ij] = Rg[ij]
        img[3, ij] = Rb[ij]
        img[4, ij] = 1.0
    end

    imgp = process_image(img; α=α, β=β)

    return colorview(RGBA, imgp)
end






function meshcube(o=Vec3f(0), sizexyz=Vec3f(1))
    uvs = map(v -> v ./ (3, 2), Vec2f[
        (0, 0), (0, 1), (1, 1), (1, 0),
        (1, 0), (1, 1), (2, 1), (2, 0),
        (2, 0), (2, 1), (3, 1), (3, 0),
        (0, 1), (0, 2), (1, 2), (1, 1),
        (1, 1), (1, 2), (2, 2), (2, 1),
        (2, 1), (2, 2), (3, 2), (3, 1),
    ])
    m = normal_mesh(Rect3f(Vec3f(-0.5) .+ o, sizexyz))
    m = GeometryBasics.Mesh(meta(coordinates(m);
            uv=uvs, normals=normals(m)), faces(m))
end





function vis_cube(hsi; cmap=:rainbow, offset=0.0, shading=false, ibounds=(1,1), jbounds=(1,1))
    size(hsi.Reflectance)

    data = PermutedDimsArray(hsi.Reflectance, (2, 3, 1))
    if ibounds != (1,1)
        imin,imax = ibounds
        data = data[imin:imax,:,:]
    end

    if jbounds != (1,1)
        jmin,jmax=jbounds
        data = data[:,jmin:jmax,:]
    end



    # make sure to clip to reasonable values
    Rmin = minimum(data)
    Rmax = maximum(data)

    cg = cgrad(cmap, range(Rmin, stop=Rmax, length=256))

    y₋ = [cg[data[i, 1, k] + offset] for i in axes(data, 1), k in axes(data, 3)][:, end:-1:1]
    y₊ = [cg[data[i, end, k] + offset] for i in axes(data, 1), k in axes(data, 3)][:, :]
    y_∅ = [RGBA(0, 0, 0, 0) for i in axes(data, 1), k in axes(data, 3)]

    x₋ = [cg[data[1, j, k] + offset] for j in axes(data, 2), k in axes(data, 3)][end:-1:1, end:-1:1]
    x₊ = [cg[data[end, j, k] + offset] for j in axes(data, 2), k in axes(data, 3)][:, end:-1:1]
    x_∅ = [RGBA(0, 0, 0, 0) for j in axes(data, 2), k in axes(data, 3)]

    z₋ = [cg[data[i, j, 300] + offset] for i in axes(data, 1), j in axes(data, 2)][end:-1:1, :]
    z₊ = getRGB(hsi)
    if ibounds != (1,1)
        imin, imax = ibounds
        z₊ = z₊[imin:imax, :]
    end

    if jbounds != (1,1)
        jmin,jmax=jbounds
        z₊ = z₊[:, jmin:jmax]
    end

    z₊ = z₊[end:-1:1, :]

    z_∅ = [RGBA(0, 0, 0, 0) for i in axes(data, 1), j in axes(data, 2)]

    m = meshcube(Vec3f(0), Vec3f(size(data)...))


    # +z, +x, +y,
    # -x, -y, -z
    img_z₊ = [
        z₊ z_∅ z_∅
        z_∅ z_∅ z_∅
    ]

    img_x₋ = [
        x_∅ x₋ x_∅
        x_∅ x_∅ x_∅
    ]

    img_x₊ = [
        rotr90(x_∅) rotr90(x_∅) rotr90(x_∅)
        rotr90(x₊) rotr90(x_∅) rotr90(x_∅)
    ]

    img_y₊ = [
        y_∅ y_∅ y_∅
        y_∅ y₊ y_∅
    ]

    img_y₋ = [
        rotl90(y_∅) rotl90(y_∅) rotl90(y₋)
        rotl90(y_∅) rotl90(y_∅) rotl90(y_∅)
    ]


    fig = Figure(; resolution=(1200, 600))
    ax = Axis3(fig[1, 1],
        aspect=:data,
        xgridvisible=false,
        ygridvisible=false,
        perspectiveness=1,
    )
    hidedecorations!(ax)
    hidespines!(ax)
    obj = mesh!(ax, m; color=img_z₊, interpolate=false, shading=shading)

    mesh!(ax, m; color=img_x₋, interpolate=false, shading=shading)
    mesh!(ax, m; color=img_x₊, interpolate=false, shading=shading)
    mesh!(ax, m; color=img_y₋, interpolate=false, shading=shading)
    mesh!(ax, m; color=img_y₊, interpolate=false, shading=shading)

    return fig
end




end
