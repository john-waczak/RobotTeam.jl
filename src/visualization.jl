function process_image(img)
    img_out = copy(img)
    # 1. normalize each band
    for i ∈ axes(img,1)
        #brighten
        α=10.0
        β=0

        img_out[i,:,:] .= clamp.(α.*img_out[i,:,:] .+ β, 0, 1)

        # normalize
        # bmin, bmax = extrema(img_out[i,:,:])
        # img_out[i,:,:] .= (img_out[i,:,:] .- bmin)./(bmax - bmin)

    end
    return img_out
end

