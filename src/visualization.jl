module Visualization

using ..ENVI: HyperspectralImage

using Images
using HDF5
using GLMakie
using GeometryBasics
using ProgressMeter
using Statistics

export process_image
export getRGB
export vis_cube
export vis_rectified_cube

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

    img = hsi.Datacube[[idx_r, idx_g, idx_b], :, :]

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

    Rr = h5["data-Δx_$(Δx)/Data"][idx_r, :, :]
    Rg = h5["data-Δx_$(Δx)/Data"][idx_g, :, :]
    Rb = h5["data-Δx_$(Δx)/Data"][idx_b, :, :]

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





function vis_cube(hsi; cmap=:jet, offset=0.1, shading=false, ibounds=(1,1), jbounds=(1,1), resolution=(800,600), azimuth=-π/4, elevation=3π/16)
    size(hsi.Reflectance)

    data = log10.(permutedims(hsi.Reflectance .+ offset, (2, 3, 1)))
    dmin, dmax = extrema(data)

    data .= (data .- dmin) ./ (dmax - dmin)

    println(extrema(data))

    if ibounds != (1, 1)
        imin, imax = ibounds
        data = data[imin:imax, :, :]
    end

    if jbounds != (1, 1)
        jmin, jmax = jbounds
        data = data[:, jmin:jmax, :]
    end



    cg = cgrad(cmap)

    y₋ = [cg[data[i, 1, k]] for i in axes(data, 1), k in axes(data, 3)][:, end:-1:1]
    y₊ = [cg[data[i, end, k]] for i in axes(data, 1), k in axes(data, 3)][:, :]
    y_∅ = [RGBA(0, 0, 0, 0) for i in axes(data, 1), k in axes(data, 3)]

    x₋ = [cg[data[1, j, k]] for j in axes(data, 2), k in axes(data, 3)][end:-1:1, end:-1:1]
    x₊ = [cg[data[end, j, k]] for j in axes(data, 2), k in axes(data, 3)][:, end:-1:1]
    x_∅ = [RGBA(0, 0, 0, 0) for j in axes(data, 2), k in axes(data, 3)]

    z₋ = [cg[data[i, j, 300]] for i in axes(data, 1), j in axes(data, 2)][end:-1:1, :]
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


    fig = Figure(; resolution=resolution, backgroundcolor=RGBA(0,0,0,0))
    ax = Axis3(
        fig[1, 1],
        aspect=:data,
        xgridvisible=false,
        ygridvisible=false,
        perspectiveness=1,
        elevation = elevation,
        azimuth = azimuth,
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



function vis_rectified_cube(Data, xs, ys, IsInbounds, λs, Δx; ibounds=(1, 1), jbounds=(1, 1), offset = 0.1, colormap=:jet, resolution=(800, 600), azimuth=3π / 4, elevation=3π / 16, colorbar=false)

    nλs = length(λs)

    println("\tGenerating RGB view...")
    ij_pixels = findall(IsInbounds)
    Ref_img = getRGB(Data, λs, ij_pixels)

    println("\tReshaping data for plotting...")
    data = PermutedDimsArray(Data[1:nλs, :, :], (2, 3, 1))

    println("\tApplying log10...")
    data = log10.(data .+ offset)

    idx_not_nan_or_inf = findall(.!(isnan.(data)) .&& .!(isinf.(data)))
    Rmin = quantile(data[idx_not_nan_or_inf], 0.1)
    Rmax = quantile(data[idx_not_nan_or_inf], 0.99)

    println("\t 0.10 quantile log10 reflectance: ", Rmin)
    println("\t 0.99 quantile log10 reflectance: ", Rmax)

    if ibounds != (1, 1)
        imin, imax = ibounds
        data = data[imin:imax, :, :]
        Ref_img = Ref_img[imin:imax, :]
    end

    if jbounds != (1, 1)
        jmin, jmax = jbounds
        data = data[:, jmin:jmax, :]
        Ref_img = Ref_img[:, jmin:jmax]
    end


    println("\tGenerating visualition...")
    fig = Figure(; resolution=resolution)
    ax = Axis3(
        fig[1, 1],
        perspectiveness=0.5,
        elevation=elevation,
        azimuth=azimuth,
        tellheight=true,
        tellwidth=true,
        # aspect=:data,
    )

    hidedecorations!(ax)
    hidespines!(ax)

    @showprogress for k in 1:length(λs)
        mr = Rect3f(Vec3f(0,0,(k-1)*Δx), Vec3f(length(xs)*Δx, length(ys)*Δx, Δx))
        mesh!(
            ax,
            mr;
            color = data[:,:,k],
            interpolate=false,
            colormap = colormap,
            colorrange= (Rmin, Rmax),
            shading=false,
            )

    end

    mr = Rect3f(Vec3f(0,0,462*Δx), Vec3f(length(xs)*Δx, length(ys)*Δx, Δx))
    mesh!(
        ax,
        mr;
        color= Ref_img,
        interpolate=false,
        shading=false,
    )

    if colorbar
        Colorbar(fig[1, 2], limits=(Rmin, Rmax), colormap=colormap, label="log10(Reflectance)", height=Relative(0.5), tellwidth=true)
    end

    return fig
end



end
