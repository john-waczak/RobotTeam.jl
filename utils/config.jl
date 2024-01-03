# basepath = "/media/jwaczak/Data/robotteam-data/raw"
# basepath = "/media/jwaczak/LabData/RobotTeam/raw/hsi"
basepath = "/Volumes/LabData/RobotTeam/raw/hsi"
@assert ispath(basepath)


# outpath = "/media/jwaczak/Data/robotteam-data/h5"
# outpath = "/media/jwaczak/LabData/RobotTeam/processed/hsi"
outpath = "/Users/johnwaczak/data/robot-team/processed/hsi"
if !ispath(outpath)
    mkpath(outpath)
end


CollectionsDict = Dict(
    "11-23" => Dict(
        "Scotty_1" => [],
        "Scotty_2" => [],
        "Scotty_3" => [],
        "Scotty_4" => [],
        "Scotty_5" => [],
    ),
    "12-09" => Dict(
        "NoDye_1" => [],
        "NoDye_2" => [],
        "Dye_1" => [],
        "Dye_2" => [],
    ),
    "12-10" => Dict(
        "NoDye_1" => [],
        "NoDye_2" => [],
        "Dye_1" => [],
        "Dye_2" => [],
    ),
    # "03-24" => Dict(
    #     "Demonstration" => [],
    #     "Demonstration_long" => [],
    # )
)


# look up bil files
for (day, runs) ∈ CollectionsDict
    for run ∈ keys(runs)
        if !ispath(joinpath(outpath, day, run))
            mkpath(joinpath(outpath, day, run))
        end

        CollectionsDict[day][run] = get_raw_file_list.(get_bil_files(joinpath(basepath, day), run))
    end
end


θ_view=30.8              # viewing angle
z_ground=292.0           # ground height (m)
isflipped=false          # flip pixel orientation
Δx = 0.10                # resampling resolution (m)
is_spec_chunked=false    # chunk HDF5 by pixel
is_band_chunked=false    # chunk HDF5 by band


# define bounding box for Scotty's Ranch
w= -97.717472
n= 33.703572
s= 33.700797
e= -97.712413




features_dict = Dict(
    :varnames => [
        ["R_"*lpad(idx, 3,"0") for idx ∈ 1:462]...,
        "roll",
        "pitch",
        "heading",
        "view_angle",
        "solar_azimuth",
        "solar_elevation",
        "solar_zenith",
        "DVI",
        "EVI",
        "GEMI",
        "GARI",
        "GCI",
        "GDVI",
        "GLI",
        "GNDVI",
        "GOSAVI",
        "GRVI",
        "GSAVI",
        "IPVI",
        "LAI",
        "LCI",
        "MNLI",
        "MSAVI2",
        "Modified Simple Ratio",
        "NLI",
        "NDVI",
        "NPCI",
        "OSAVI",
        "RDVI",
        "SAVI",
        "SR",
        "TDVI",
        "TGI",
        "VARI",
        "WDRVI",
        "AVRI",
        "MCARI",
        "MCARI2",
        "MRENDVI",
        "MRESER",
        "MTVI",
        "RENDVI",
        "TCARI",
        "TVI",
        "VREI1",
        "VREI2",
        "VREI3",
        "PRI",
        "SIPI",
        "SIPI1",
        "PSRI",
        "ARI1",
        "ARI2",
        "CRI1",
        "CRI2",
        "NDWI1",
        "NDWI2",
        "MNDWI",
        "WBI",
        "ACI",
        "MARI",
        "MSI",
        "MTCI",
        "NDII",
        "NDRE",
        "RGRI",
        "RVSI",
        "yaw_minus_azimuth",
        "Σrad",
        "Σdownwelling",
    ],
    :printnames => [
        ["Reflectance Band "*lpad(idx, 3,"0") for idx ∈ 1:462]...,
        "Roll",
        "Pitch",
        "Heading",
        "Viewing Angle",
        "Solar Azimuth",
        "Solar Elevation",
        "Solar Zenith",
        "Difference Vegetation Index",
        "Enhanced Vegetation Index",
        "Global Environmental Monitoring Index",
        "Green Atmospherically Resistant Index",
        "Green Chlorophyll Index",
        "Green Difference Vegetation Index",
        "Green Leaf Index",
        "Green Normalized Difference Vegetation Index",
        "Green Optimized Soil Adjusted Vegetation Index",
        "Green Ratio Vegetation Index",
        "Green Soil Adjusted Vegetation Index",
        "Infrared Percentage Vegetation Index",
        "Leaf Area Index",
        "Leaf Chlorphyll Index",
        "Modified Non-Linear Index",
        "Modified Soil Adjusted Vegetation Index 2",
        "Modified Simple Ratio",
        "Non-Linear Index",
        "Normalized Difference Vegetation Index",
        "Normalized Pigment Chlorophyll Index",
        "Optimized Soil Ajusted Vegetation Index",
        "Renormalized Difference Vegetation Index",
        "Simple Ratio",
        "Transformed Difference Vegetation Index",
        "Triangular Greenness Index",
        "Visible Atmospherically Resistant Index",
        "Wide Dynamic Range Vegetation Index",
        "Atmospherically Resistant Vegetation Index",
        "Modified Chlorophyll Absorption Ratio Index",
        "Modified Chlorophyll Absorption Ratio Index Improved",
        "Modified Red Edge Normalized Difference Vegetation Index",
        "Modified Red Edge Simple Ratio",
        "Modified Triangular Vegetation Index",
        "Red Edge Normalized Difference Vegetation Index",
        "Transformed Chlorophyll Absorption Reflectance Index",
        "Triangular Vegetation Index",
        "Vogelmann Red Edge Index 1",
        "Vogelmann Red Edge Index 2",
        "Vogelmann Red Edge Index 3",
        "Photochemical Reflectance Index",
        "Structure Insensitive Pigment Index",
        "Structure Independent Pigment Index",
        "Plant Senescence Reflectance Index",
        "Anthocyanin Reflectance Index 1",
        "Anthocyanin Reflectance Index 2",
        "Carotenoid Reflectance Index 1",
        "Carotenoid Reflectance Index 2",
        "Normalized Difference Water Index 1",
        "Normalized Difference Water Index 2",
        "Modified Normalized Difference Water Index",
        "Water Band Index",
        "Anthocyanin Content Index",
        "Modified Anthocyanin Reflectance Index",
        "Moisture Stress Index",
        "MERIS Terrestrial Chlorophyll Index",
        "Normalized Difference Infrared Index",
        "Normalized Difference Red Edge",
        "Red Green Ratio Index",
        "Red Edge Vegetation Stress Index",
        "Heading - Solar Azimuth",
        "Total Pixel Intensity",
        "Total Downwelling Intensity",
    ],
)
