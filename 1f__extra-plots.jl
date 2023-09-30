using CairoMakie
using MintsMakieRecipes
using XLSX
using CSV, DataFrames
using DataInterpolations

set_theme!(mints_theme)
update_theme!(
    figure_padding=30,
    Axis = (
        xticklabelsize=20,
        yticklabelsize=20,
        xlabelsize=22,
        ylabelsize=22,
        titlesize=25,
    ),
    Colorbar = (
        ticklabelsize=20,
        labelsize=22
    )
)




# spectral data from the USGS
# https://landsat.usgs.gov/spectral-characteristics-viewer

basepath = "./assets/satellite-passbands"
@assert ispath(basepath)

readdir(basepath)

fnames = [
    "L8_OLI_RSR.xlsx",
    "MODIS_PFM_IB_OOB_RSR_merged.xlsx",
    "S2-SRF_COPE-GSEG-EOPG-TN-15-0007_3.1.xlsx",
]

files = joinpath.(basepath, fnames)

@assert all(isfile.(files))





# --------------------------------------------------------
# Sentinel 2
# --------------------------------------------------------

out_name = ["λ", "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A", "B9", "B10", "B11", "B12"]

df_s2a = DataFrame(XLSX.readtable(files[end], "Spectral Responses (S2A)"))
df_s2b = DataFrame(XLSX.readtable(files[end], "Spectral Responses (S2B)"))

for i ∈ 1:ncol(df_s2a)
    rename!(df_s2a, i => out_name[i])
    rename!(df_s2b, i => out_name[i])

    if i == 1
        df_s2a[!,i] = Int.(df_s2a[!,i])
        df_s2b[!,i] = Int.(df_s2b[!,i])
    else
        df_s2a[!,i] = Float64.(df_s2a[!,i])
        df_s2b[!,i] = Float64.(df_s2b[!,i])
    end
end


# --------------------------------------------------------
# Landsat 8
# --------------------------------------------------------

λs = df_s2a.λ

df_ls8 = DataFrame(
    :λ => λs
)

xf = XLSX.readxlsx(files[1])
XLSX.sheetnames(xf)

let
    df_b = DataFrame(XLSX.readtable(files[1], "CoastalAerosol"));

    rsr_out = zeros(length(λs));

    for row ∈ eachrow(df_b)
        idx=findfirst(row.Wavelength .== λs);
        rsr_out[idx] = row["Band 1"];
    end

    df_ls8[!, "B1"] = rsr_out;
end


let
    df_b = DataFrame(XLSX.readtable(files[1], "Blue"));

    rsr_out = zeros(length(λs));

    for row ∈ eachrow(df_b)
        idx=findfirst(row.Wavelength .== λs);
        rsr_out[idx] = row["Band 2"];
    end

    df_ls8[!, "B2"] = rsr_out;
end

let
    df_b = DataFrame(XLSX.readtable(files[1], "Green"));

    rsr_out = zeros(length(λs));

    for row ∈ eachrow(df_b)
        idx=findfirst(row.Wavelength .== λs);
        rsr_out[idx] = row["Band 3"];
    end

    df_ls8[!, "B3"] = rsr_out;
end

let
    df_b = DataFrame(XLSX.readtable(files[1], "Red"));

    rsr_out = zeros(length(λs));

    for row ∈ eachrow(df_b)
        idx=findfirst(row.Wavelength .== λs);
        rsr_out[idx] = row["Band 4"];
    end

    df_ls8[!, "B4"] = rsr_out;
end

let
    df_b = DataFrame(XLSX.readtable(files[1], "NIR"));

    rsr_out = zeros(length(λs));

    for row ∈ eachrow(df_b)
        idx=findfirst(row.Wavelength .== λs);
        rsr_out[idx] = row["Band 5"];
    end

    df_ls8[!, "B5"] = rsr_out;
end

let
    df_b = DataFrame(XLSX.readtable(files[1], "Cirrus"));

    rsr_out = zeros(length(λs));

    for row ∈ eachrow(df_b)
        idx=findfirst(row.Wavelength .== λs);
        rsr_out[idx] = row["Band 9"];
    end

    df_ls8[!, "B9"] = rsr_out;
end

let
    df_b = DataFrame(XLSX.readtable(files[1], "SWIR1"));

    rsr_out = zeros(length(λs));

    for row ∈ eachrow(df_b)
        idx=findfirst(row.Wavelength .== λs);
        rsr_out[idx] = row["Band 6"];
    end

    df_ls8[!, "B6"] = rsr_out;
end

let
    df_b = DataFrame(XLSX.readtable(files[1], "SWIR2"));

    rsr_out = zeros(length(λs));

    for row ∈ eachrow(df_b)
        idx=findfirst(row.Wavelength .== λs);
        rsr_out[idx] = row["Band 7"];
    end

    df_ls8[!, "B7"] = rsr_out;
end

let
    df_b = DataFrame(XLSX.readtable(files[1], "Pan"));

    rsr_out = zeros(length(λs));

    for row ∈ eachrow(df_b)
        idx=findfirst(row.Wavelength .== λs);
        rsr_out[idx] = row["Band 8"];
    end

    df_ls8[!, "B8"] = rsr_out;
end


# --------------------------------------------------------
# MODIS
# --------------------------------------------------------
df_modis = DataFrame()
df_modis.λ = λs

xf = XLSX.readxlsx(files[2])
nm = XLSX.sheetnames(xf)

data = XLSX.readdata(files[2], nm[1]*"!A3:BT320")

i_band = 1
for j∈1:2:size(data,2)
    D = data[:,j:j+1]
    D = D[.!(ismissing.(D[:,1])), :]
    D = Float64.(D)

    λs_modis = 1000.0 * D[:,1]
    rsr = D[:,2]

    itp = LinearInterpolation(rsr, λs_modis)
    rsr_out = itp.(λs)
    rsr_out[findall(λs .< λs_modis[1])] .= 0.0
    rsr_out[findall(λs .> λs_modis[end])] .= 0.0

    df_modis[!, "B$(i_band)"] = rsr_out
    i_band += 1
end

# turn this into a dictionary with nonzero wavelengths only
size_in_inches = (16/4, 9/4)
dpi = 300
size_in_pixels = size_in_inches .* dpi

fig = Figure(resolution=size_in_pixels);
g = fig[1,1] = GridLayout()

ax11 = Axis(g[1,1], ylabel="Sentinel 2A",
            yticklabelsvisible=false, yticksvisible=false,
            xticksvisible=false, xticklabelsvisible=false
            );

ax12 = Axis(g[1,2],
            yticklabelsvisible=false, yticksvisible=false,
            xticksvisible=false, xticklabelsvisible=false
            );

ax21 = Axis(g[2,1], ylabel="Sentinel 2B",
            yticklabelsvisible=false, yticksvisible=false,
            xticksvisible=false, xticklabelsvisible=false
            );

ax22 = Axis(g[2,2],
            yticklabelsvisible=false, yticksvisible=false,
            xticksvisible=false, xticklabelsvisible=false
            );

ax31 = Axis(g[3,1], ylabel="Landsat 8",
            yticklabelsvisible=false, yticksvisible=false,
            xticksvisible=false, xticklabelsvisible=false
            );

ax32 = Axis(g[3,2],
            yticklabelsvisible=false, yticksvisible=false,
            xticksvisible=false, xticklabelsvisible=false
            );

ax41 = Axis(g[4,1], ylabel="Modis",
            xlabel="VNIR Wavelengths (nm)",
            yticklabelsvisible=false, yticksvisible=false,
            );

ax42 = Axis(g[4,2],
            xlabel="SWIR Wavelengths (nm)",
            yticksvisible=false, yticklabelsvisible=false
            );


ylims!(ax11, 0.0, 1.0)
ylims!(ax21, 0.0, 1.0)
ylims!(ax31, 0.0, 1.0)
ylims!(ax41, 0.0, 1.0)

ylims!(ax12, 0.0, 1.0)
ylims!(ax22, 0.0, 1.0)
ylims!(ax32, 0.0, 1.0)
ylims!(ax42, 0.0, 1.0)



fig

linkyaxes!(ax11, ax12)
linkyaxes!(ax21, ax22)
linkyaxes!(ax31, ax32)
linkyaxes!(ax41, ax42)

linkxaxes!(ax11, ax21)
linkxaxes!(ax21, ax31)
linkxaxes!(ax31, ax41)

linkxaxes!(ax12, ax22)
linkxaxes!(ax22, ax32)
linkxaxes!(ax32, ax42)


fig


# Plot Sentinel 2A
for i ∈ 2:ncol(df_s2a)
    idx_VNIR = findall(df_s2a.λ .≤ 1000.0)
    idx_SWIR = findall(df_s2a.λ .> 1000.0)

    band!(ax11, df_s2a.λ[idx_VNIR], df_s2a[idx_VNIR,i], zeros(length(idx_VNIR)), color=(mints_colors[1], 0.2))
    lines!(ax11, df_s2a.λ[idx_VNIR], df_s2a[idx_VNIR,i], color=(mints_colors[1], 0.9))

    band!(ax12, df_s2a.λ[idx_SWIR], df_s2a[idx_SWIR,i], zeros(length(idx_SWIR)), color=(mints_colors[1], 0.2))
    lines!(ax12, df_s2a.λ[idx_SWIR], df_s2a[idx_SWIR,i], color=(mints_colors[1], 0.9))

end

for i ∈ 2:ncol(df_s2a)
    idx_VNIR = findall(df_s2b.λ .≤ 1000.0)
    idx_SWIR = findall(df_s2b.λ .> 1000.0)

    band!(ax21, df_s2b.λ[idx_VNIR], df_s2b[idx_VNIR,i], zeros(length(idx_VNIR)), color=(mints_colors[2], 0.2))
    lines!(ax21, df_s2b.λ[idx_VNIR], df_s2b[idx_VNIR,i], color=(mints_colors[2], 0.9))

    band!(ax22, df_s2b.λ[idx_SWIR], df_s2b[idx_SWIR,i], zeros(length(idx_SWIR)), color=(mints_colors[2], 0.2))
    lines!(ax22, df_s2b.λ[idx_SWIR], df_s2b[idx_SWIR,i], color=(mints_colors[2], 0.9))
end

for i ∈ 2:ncol(df_ls8)
    idx_VNIR = findall(df_ls8.λ .≤ 1000.0)
    idx_SWIR = findall(df_ls8.λ .> 1000.0)

    band!(ax31, df_ls8.λ[idx_VNIR], df_ls8[idx_VNIR,i], zeros(length(idx_VNIR)), color=(mints_colors[3], 0.2))
    lines!(ax31, df_ls8.λ[idx_VNIR], df_ls8[idx_VNIR,i], color=(mints_colors[3], 0.9))

    band!(ax32, df_ls8.λ[idx_SWIR], df_ls8[idx_SWIR,i], zeros(length(idx_SWIR)), color=(mints_colors[3], 0.2))
    lines!(ax32, df_ls8.λ[idx_SWIR], df_ls8[idx_SWIR,i], color=(mints_colors[3], 0.9))
end

for i ∈ 2:ncol(df_modis)
    idx_VNIR = findall(df_modis.λ .≤ 1000.0)
    idx_SWIR = findall(df_modis.λ .> 1000.0)

    band!(ax41, df_modis.λ[idx_VNIR], df_modis[idx_VNIR,i], zeros(length(idx_VNIR)), color=(mints_colors[4], 0.2))
    lines!(ax41, df_modis.λ[idx_VNIR], df_modis[idx_VNIR,i], color=(mints_colors[4], 0.9))

    band!(ax42, df_modis.λ[idx_SWIR], df_modis[idx_SWIR,i], zeros(length(idx_SWIR)), color=(mints_colors[4], 0.2))
    lines!(ax42, df_modis.λ[idx_SWIR], df_modis[idx_SWIR,i], color=(mints_colors[4], 0.9))
end


Label(g[:,:,Top()], "Relative Spectral Response", fontsize=25, font=:bold, padding=(0,5,10,0))

colsize!(g, 1, Relative(7/12))


fig



save("./paper/assets/passbands.svg", fig, pt_per_unit=1)
save("./paper/assets/passbands.png", fig, pt_per_unit=1)
save("./paper/assets/passbands.eps", fig, pt_per_unit=1)
save("./paper/assets/passbands.pdf", fig, pt_per_unit=1)



