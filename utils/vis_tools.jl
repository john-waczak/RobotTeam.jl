using CondaPkg
CondaPkg.add("contextily")
using PythonCall
using Images

export get_background_satmap


ctx = pyimport("contextily")

struct SatMap
    w
    e
    s
    n
    img
end


"""
    get_background_satmap(w::Float64, e::Float64, s::Float64, n::Float64; out_name::String="Scotty)


Grab Esri World Imagery tiles for region with corners (w,n), (e,s) in longitude and latitude.
Saves resulting image to `outpath`

**Note:** result images are saved in Web-Mercator projection by default. See `WebMercatorfromLLA` and `LLAfromWebMercator` from `Geodesy.jl` for conversion details.
"""
function get_background_satmap(w::Float64, e::Float64, s::Float64, n::Float64; out_name::String="Scotty")
    # ctx = pyimport("contextily")
    tiff, ext = ctx.bounds2raster(w, s, e, n, out_name*".tiff", source=ctx.providers["Esri"]["WorldImagery"], ll=true)
    warped_tiff, warped_ext = ctx.warp_tiles(tiff, ext, "EPSG:4326")

    warped_ext = pyconvert(Vector{Float64}, warped_ext)
    tiff_array = permutedims(pyconvert(Array{Int}, warped_tiff)./255, (3,1,2))
    tiff_img = colorview(RGBA, tiff_array)

    tiff_img = rotr90(tiff_img)

    return SatMap(
        warped_ext[1],
        warped_ext[2],
        warped_ext[3],
        warped_ext[4],
        tiff_img
    )
end



# metric bounds
metric_bounds = Dict(
    "mNDWI" => (-1, 1),
    "NDVI" => (-1 ,1),
    "SR" => (0, 30),
    "EVI" => (-1, 1),
    "AVRI" => (-1, 1),
    "NDVI_705" => (-1,1),
    "MSR_705" => (0, 30),
    "MNDVI" => (-1, 1),
    "VOG1" => (0, 20),
    "VOG2" => (0, 20),
    "VOG3" => (0, 20),
    "PRI" => (-1,1),
    "SIPI" => (0, 2),
    "PSRI" => (-1, 1),
    "CRI1" => (0, 15),
    "CRI2" => (0, 15),
    "ARI1" => (0, 0.2),
    "ARI2" => (0, 0.2),
    "WBI" => (0.5, 1.5),
    "MCRI" => (0, 15),
    "TCARI" => (0, 15)
)
