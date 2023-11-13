module Visualization

using ..ENVI: HyperspectralImage

using Images
using HDF5
# using GLMakie
using ProgressMeter
using Statistics

export process_image
export getRGB
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
