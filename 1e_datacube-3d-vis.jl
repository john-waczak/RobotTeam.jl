using RobotTeam
using HDF5
using BenchmarkTools
using GLMakie, GeometryBasics
using MintsMakieRecipes

set_theme!(mints_theme)

include("utils/vis_tools.jl")
include("utils/config.jl")

basepath

files = get_raw_file_list.(get_bil_files(basepath, "Scotty_1"))

f = files[1]


hsi = HyperspectralImage(f.bilpath, f.bilhdr, f.lcfpath, f.timespath, f.specpath, f.spechdr; isflipped=true);


size(hsi.Reflectance)

data = PermutedDimsArray(hsi.Reflectance, (2,3,1))

# make sure to clip to reasonable values
Rmin = minimum(data)
Rmax = maximum(data)
cg = cgrad(:inferno, range(Rmin, stop=Rmax, length=256))


size(data)

heatmap(z₊)

heatmap(data[1,:,:])

size(y₋)
size(x₋)

y₋ = [cg[data[i,1,k]] for i in axes(data,1), k in axes(data,3)][:, end:-1:1]
y₊ = [cg[data[i,end,k]] for i in axes(data,1), k in axes(data,3)][:, :]
y_∅ = [RGBA(0,0,0,0) for i in axes(data,1), k in axes(data,3)]

x₋ = [cg[data[1,j,k]] for j in axes(data,2), k in axes(data,3)][:,end:-1:1]
x₊ = [cg[data[end,j,k]] for j in axes(data,2), k in axes(data,3)][:,end:-1:1]
x_∅ = [RGBA(0,0,0,0) for j in axes(data,2), k in axes(data,3)]

z₋ = [cg[data[i,j,300]] for i in axes(data,1), j in axes(data,2)][end:-1:1,:]
z₊ = getRGB(hsi)[end:-1:1,:]
z_∅ = [RGBA(0,0,0,0) for i in axes(data,1), j in axes(data,2)]

function meshcube(o=Vec3f(0), sizexyz = Vec3f(1))
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
                                 uv = uvs, normals = normals(m)), faces(m))
end

m = meshcube(Vec3f(0), Vec3f(size(data)...))


# +z, +x, +y,
# -x, -y, -z

img_z₊ = [
    z₊ z_∅ z_∅;
    z_∅ z_∅ z_∅;
]

img_x₋ = [
    x_∅ x₋ x_∅;
    x_∅ x_∅ x_∅;
]

img_x₊ = [
    rotr90(x_∅) rotr90(x_∅) rotr90(x_∅);
    rotr90(x₊) rotr90(x_∅) rotr90(x_∅);
]

img_y₊ = [
    y_∅ y_∅ y_∅;
    y_∅ y₊ y_∅;
]

img_y₋ = [
    rotl90(y_∅) rotl90(y_∅) rotl90(y₋);
    rotl90(y_∅) rotl90(y_∅) rotl90(y_∅);
]



fig = Figure(; resolution=(1200, 600))
ax = Axis3(fig[1,1],
           aspect=:data,
           xgridvisible=false,
           ygridvisible=false,
           perspectiveness=1,
           )
hidedecorations!(ax)
obj = mesh!(ax, m; color = img_z₊, interpolate=false)

mesh!(ax, m; color=img_x₋, interpolate=false)
mesh!(ax, m; color=img_x₊, interpolate=false)
mesh!(ax, m; color=img_y₋, interpolate=false)
mesh!(ax, m; color=img_y₊, interpolate=false)


fig



heatmap(z₊)

# xs, ys, isnorth, zone, Longitudes, Latitudes, IsInbounds, varnames, printnames, λs, Data_μ = resample_datacube_fast(hsi; Δx=Δx)
