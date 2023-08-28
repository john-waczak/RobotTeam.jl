using RobotTeam
using HDF5
using BenchmarkTools
using GLMakie, GeometryBasics
using MintsMakieRecipes
using ProgressMeter

set_theme!(mints_theme)

include("utils/vis_tools.jl")
include("utils/config.jl")



basepath

# files = get_raw_file_list.(get_bil_files(joinpath(basepath, "11-23"), "Scotty_1"))
# f = files[1]



files = get_raw_file_list.(get_bil_files(joinpath(basepath, "12-09"), "Dye_1"))
f = files[6]


hsi = HyperspectralImage(f.bilpath, f.bilhdr, f.lcfpath, f.timespath, f.specpath, f.spechdr; isflipped=true);
img = getRGB(hsi)
fig, ax, hm = heatmap(img)

size(img)


fig = vis_cube(hsi; cmap=:batlow, offset=0.0, ibounds=(250,1600))
save("demo-cube.png", fig)



# now let's visualize the georectified cube

xs, ys, isnorth, zone, Longitudes, Latitudes, IsInbounds, varnames, printnames, λs, Data_μ = resample_datacube_fast(hsi; Δx=Δx)


data = PermutedDimsArray(Data_μ[1:462,:,:], (2, 3, 1))



fig = Figure(; resolution=(1200,600))
ax = Axis3(fig[1, 1],
           # aspect=:data,
           #xgridvisible=false,
           #ygridvisible=false,
           perspectiveness=1,
           )
hidedecorations!(ax)
hidespines!(ax)

@showprogress for k in 1:length(hsi.λs)
    mr = Rect3f(Vec3f(0,0,(k-1)*Δx), Vec3f(length(xs)*Δx, length(ys)*Δx, Δx))
    mesh!(
        ax,
        mr;
        color = data[:,:,k],
        interpolate=false,
        colormap = :rainbow,
        )

end

ij_pixels = findall(IsInbounds)
Ref_img = getRGB(Data_μ, λs, ij_pixels)

mr = Rect3f(Vec3f(0,0,462*Δx), Vec3f(length(xs)*Δx, length(ys)*Δx, Δx))
mesh!(
    ax,
    mr;
    color= Ref_img,
    interpolate=false,
    colormap = :rainbow,
)

save("demo-rectified-cube.png", fig)

fig


