module Visualization

using ..ENVI: HyperspectralImage

using Images
using HDF5


export process_image
export generateRGB!
export getRGB

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







end
