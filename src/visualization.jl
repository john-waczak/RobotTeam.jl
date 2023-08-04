module Visualization

using Images
using HDF5


export process_image
export generateRGB!
export getRGB

function process_image(img;α=10.0,β=0.0)
    # see example: https://www.satmapper.hu/en/rgb-images/
    img_out = copy(img)
    # 1. normalize each band
    for i ∈ axes(img,1)
        #brighten
        img_out[i,:,:] .= clamp.(α.*img_out[i,:,:] .+ β, 0, 1)

        # normalize
        # bmin, bmax = extrema(img_out[i,:,:])
        # img_out[i,:,:] .= (img_out[i,:,:] .- bmin)./(bmax - bmin)

    end
    return img_out
end


"""
    generateRGB!(h5::HDF5.File)

Given a HSI stored in HDF5 file `h5`, generate an RGB representation of the reflectance data by selecting specific wavelengths for the red, green, and blue bins.
"""
function generateRGB!(h5::HDF5.File; λred=630.0, λgreen=532.0, λblue=465.0)
    println("\tfetching wavelengths")
    λs = read(h5["raw/reflectance/wavelengths"])
    idx_r = argmin(abs.(λs .- λred))
    idx_g = argmin(abs.(λs .- λgreen))
    idx_b = argmin(abs.(λs .- λblue))

    println("\tfetching reflectances")
    img = permutedims(read(h5["raw/reflectance/reflectance"])[[idx_r, idx_g, idx_b],:,:], (1,3,2))

    println("\tparsing as image")
    imgp = process_image(img)

    g = create_group(h5["raw"], "RGB")
    g["RGB"] = imgp
end


function getRGB(h5::HDF5.File)
    return colorview(RGB, read(h5["raw/RGB/RGB"]))
end





end
