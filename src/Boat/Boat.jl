module Boat

using ..Georectification: bump_to_nearest_Δx

using CSV, DataFrames
using Dates
using ProgressMeter
using ProgressMeter
using DataInterpolations
using Geodesy


# we will need to parse NMEA encoded gps/sensor strings
# see this site: http://aprs.gids.nl/nmea/


export importAirMar, importCOM1, importCOM2, importCOM3, importLISST, importNMEA, processBoatFiles, combine_boat_dfs




"""
    importAIRMar(path::String)

Read data from AIRMar style text files and return DataFrames with their data.
"""
function importAirMar(path::String)
    GPGGA = Dict("datetime"=>DateTime[],
                 "utc_dt" =>DateTime[],
                 "unix_dt" =>Float64[],
                 "latitude"=>Float64[],
                 "longitude"=>Float64[],
                 "altitude"=>Float64[],
                 )

    GPVTG = Dict("utc_dt" =>DateTime[],
                 "unix_dt"=>Float64[],
                 "speedInKmph"=>Float64[],
                 "speedInKnots"=>Float64[],
                 )

    #     WIMDA = Dict()
    #     TIROT = Dict()
    #     YXXDR = Dict()
    #     GPZDA = Dict()
    #     WIMWV = Dict()
    #     YXXDR = Dict()

    for line ∈ readlines(path)
        splitline = split(line, ",")
        date = splitline[1]
        time = splitline[2]
        dt = DateTime(date*" "*time, "m/d/y H:M:S") + Year(2000)
        unix_dt = datetime2unix.(dt)

        id = splitline[3]
        if occursin("GPGGA", id) && !any(isempty(s) for s ∈ splitline[1:12])
            push!(GPGGA["datetime"], dt)

            # parse utc time
            utc_raw = splitline[4]
            utc_h = utc_raw[1:2]
            utc_m = utc_raw[3:4]
            utc_s = utc_raw[5:end]
            # make sure to double check this on old data
            # utc_dt = DateTime(date*" "*utc_h*":"*utc_m*":"*utc_s, "m/d/y H:M:S.s")
            utc_dt = DateTime(date*" "*utc_h*":"*utc_m*":"*utc_s, "m/d/y H:M:S.s") + Year(2000)
            push!(GPGGA["utc_dt"], utc_dt)
            push!(GPGGA["unix_dt"], unix_dt)
            # parse longitude. NOTE: be sure to keep full precision
            lat = split(splitline[5], ".")
            deg = parse(Float64, lat[1][1:end-2])
            min = parse(Float64, lat[1][end-1:end]*"."*lat[2])
            lat = deg + min/60

            # check for north vs south
            if splitline[6] == "S"
                lat = -lat
            else
                lat = lat
            end
            push!(GPGGA["latitude"], lat)

            # parse latitude. NOTE: be sure to keep full precision
            long = split(splitline[7], ".")
            deg = parse(Float64, long[1][1:end-2])
            min = parse(Float64, long[1][end-1:end]*"."*long[2])
            long = deg + min/60

            if splitline[8] == "W"
                long = -long
            else
                long = long
            end
            push!(GPGGA["longitude"], long)


            # parse altitude.
            alt = parse(Float64, splitline[12])
            push!(GPGGA["altitude"], alt)
        end


        if occursin("GPVTG", id) && !any(isempty(s) for s ∈ splitline[8:11])
            kmh = parse(Float64, splitline[10])
            knots = parse(Float64, splitline[8])

            push!(GPVTG["speedInKmph"], kmh)
            push!(GPVTG["speedInKnots"], knots)
            push!(GPVTG["utc_dt"], dt)
            push!(GPVTG["unix_dt"], unix_dt)
        end


    end
    return DataFrame(GPGGA), DataFrame(GPVTG)
end


"""
    importNMEA(path::String)

Read an NMEA text file to a DataFrame.

# Fields
- **utc_dt**: the utc timestamp for each datum
- **unix_dt**: unix epoch time corresponding to reported utc time
- **latitude**: Negative for Southern hemisphere
- **longitude**: Negative for Western hemisphere
"""
function importNMEA(path::String)
    returnDict = Dict("utc_dt" =>DateTime[],
                      "unix_dt"=>Float64[],
                      "latitude"=>Float64[],
                      "longitude"=>Float64[],
                      )

    for line in readlines(path)
        if occursin("GNRMC", line) || occursin("GPRMC", line)
            splitline = split(line, ",")
            if !any(isempty(s) for s ∈ splitline[1:3]) && !any(isempty(s) for s ∈ splitline[6:9])

                date = splitline[1]
                time = splitline[2]
                dt = DateTime(date*" "*time, "m/d/y H:M:S") + Year(2000)
                unix_dt = datetime2unix(dt)

                # date = splitline[1]
                # time = splitline[4]
                # dt = DateTime(date*" "*time, "m/d/y HHMMSS.ss") + Year(2000)
                # unix_dt = datetime2unix(dt)

                # parse longitude. NOTE: be sure to keep full precision
                # note: format is dddmm.mmmm. See https://stackoverflow.com/questions/6619377/how-to-get-whole-and-decimal-part-of-a-number
                lat = split(splitline[6], ".")
                deg = parse(Float64, lat[1][1:end-2])
                min = parse(Float64, lat[1][end-1:end]*"."*lat[2])
                lat = deg + min/60

                # check for north vs south
                if splitline[7] == "S"
                    lat = -lat
                else
                    lat = lat
                end

                # parse latitude. NOTE: be sure to keep full precision
                long = split(splitline[8], ".")
                deg = parse(Float64, long[1][1:end-2])
                min = parse(Float64, long[1][end-1:end]*"."*long[2])
                long = deg + min/60

                if splitline[9] == "W"
                    long = -long
                else
                    long = long
                end

                push!(returnDict["utc_dt"], dt)
                push!(returnDict["unix_dt"], unix_dt)
                push!(returnDict["latitude"],lat)
                push!(returnDict["longitude"], long)

            end
        end
    end
    return DataFrame(returnDict)
end



"""
    importCOM1(path::String)

Read an COM1 textfile into a DataFrame.

# Fields
- **utc_dt**: utc timestamp [datetime]
- **unix_dt**: unix epoch time [seconds]
- **Temp3488**: Temperature [°C]
- **pH**: [dimensionless]
- **SpCond**: [μS/cm]
- **Turb3488**: turbidity [FNU]
- **Br**: Bromine [mg/l]
- **Ca**: Calcium [mg/l]
- **Cl**: Chlorine [mg/l]
- **Na**: Sodium [mg/l]
- **NO3**: Nitrate [mg/l]
- **NH4**: Ammonium [mg/l]
- **HDO**: Heavy Water [mg/l]
- **HDO_percent**: percentage of HDO/H20 [% Sat]
- **pH_mV**: pH probe raw [mV] ??? Is this right?
- **Salinity3488**: salinity in practical salinity scale [PSS]
- ** TDS**: total dissolved solids [mg/l]
"""
function importCOM1(path::String)
    COM1 = Dict("utc_dt" => DateTime[],
                "unix_dt" => Float64[],
                "Temp3488" => Float64[],
                "pH"=>Float64[],
                "SpCond" => Float64[],
                "Turb3488" => Float64[],
                "Br" => Float64[],
                "Ca" => Float64[],
                "Cl" => Float64[],
                "Na" => Float64[],
                "NO3" => Float64[],
                "NH4" => Float64[],
                "HDO" => Float64[],
                "HDO_percent" => Float64[],
                "pH_mV" => Float64[],
                "Salinity3488" => Float64[],
                "TDS" => Float64[],
                )

    for line ∈ readlines(path)
        splitline = split(line, ",")
        if length(splitline) == 20
            # push the datetime
            date = splitline[1]
            time = splitline[2]
            dt = DateTime(date*" "*time, "m/d/y H:M:S") + Year(2000)
            unix_dt = datetime2unix(dt)
            push!(COM1["utc_dt"], dt)
            push!(COM1["unix_dt"], unix_dt)

            # Temp3488
            temp = parse(Float64, splitline[6])
            push!(COM1["Temp3488"], temp)

            # pH
            pH = parse(Float64, splitline[7])
            push!(COM1["pH"], pH)

            # SpCond
            spcond = parse(Float64, splitline[8])
            push!(COM1["SpCond"], spcond)


            # Turb3488
            turb = parse(Float64, splitline[9])
            push!(COM1["Turb3488"], turb)

            # Br
            br = parse(Float64, splitline[10])
            push!(COM1["Br"], br)

            # Ca
            ca = parse(Float64, splitline[11])
            push!(COM1["Ca"], ca)

            # Cl
            cl = parse(Float64, splitline[12])
            push!(COM1["Cl"], cl)

            # Na
            na = parse(Float64, splitline[13])
            push!(COM1["Na"], na)

            #NO3
            no3 = parse(Float64, splitline[14])
            push!(COM1["NO3"], no3)

            # NH4
            nh4 = parse(Float64, splitline[15])
            push!(COM1["NH4"], nh4)

            # HDO
            hdo = parse(Float64, splitline[16])
            push!(COM1["HDO"], hdo)

            # HDO_percent
            hdo_per = parse(Float64, splitline[17])
            push!(COM1["HDO_percent"], hdo_per)

            # pH_mV
            ph_mv = parse(Float64, splitline[18])
            push!(COM1["pH_mV"], ph_mv)

            # Salinity3488
            sal = parse(Float64, splitline[19])
            push!(COM1["Salinity3488"], sal)

            # TDS
            tds = parse(Float64, splitline[20])
            push!(COM1["TDS"], tds)



        end
    end
    return DataFrame(COM1)
end


"""
    importCOM2(path::String)

Read a COM2 textfile into a DataFrame.

# Fields
- **utc_dt**: datetime
- **unix_dt**: unix epoch time [seconds]
- **Temp3489**: Temperature [°C]
- **bg**: [ppb]
- **bgm**: [ppb]
- **CDOM**: colored dissolved organic matter [ppb]
- **Chl**: [μg/l]
- **ChlRed**: [μg/l]
- **Turb3489**: turbidity [FNU]
"""
function importCOM2(path::String)
    COM2 = Dict("utc_dt"=>DateTime[],
                "unix_dt"=>Float64[],
                "Temp3489"=>Float64[],
                "bg" => Float64[],
                "bgm" => Float64[],
                "CDOM" => Float64[],
                "Chl" => Float64[],
                "ChlRed" => Float64[],
                "Turb3489" => Float64[]
                )
    for line ∈ readlines(path)
        splitline = split(line, ",")

        if length(splitline) == 12
            # push the datetime
            date = splitline[1]
            time = splitline[2]
            dt = DateTime(date*" "*time, "m/d/y H:M:S") + Year(2000)
            unix_dt = datetime2unix(dt)
            push!(COM2["utc_dt"], dt)
            push!(COM2["unix_dt"], unix_dt)

            # Temp3488
            temp = parse(Float64, splitline[6])
            push!(COM2["Temp3489"], temp)

            # bg
            bg = parse(Float64, splitline[7])
            push!(COM2["bg"], bg)

            # bgm
            bgm = parse(Float64, splitline[8])
            push!(COM2["bgm"], bgm)

            # CDOM
            cdom = parse(Float64, splitline[9])
            push!(COM2["CDOM"], cdom)

            # chl
            chl = parse(Float64, splitline[10])
            push!(COM2["Chl"], chl)

            # ChlRed
            chlred = parse(Float64, splitline[11])
            push!(COM2["ChlRed"], chlred)

            # Turb3489
            turb = parse(Float64, splitline[12])
            push!(COM2["Turb3489"], turb)
        end
    end
    return DataFrame(COM2)
end



"""
    importCOM3(path::String)

Read a COM3 text file into a DataFrame

# Fields
- **utc_dt**: datetime
- **unix_dt**: unix epoch time [seconds]
- **Temp3490**: temperature [°C]
- **CO**: Crude Oil [ppb]
- **OB**: [ppb]
- **RefFuel**: [ppb]
- **TRYP**: [ppb]
- **Turb3490**: Turbidity [FNU]
- **Salinity3490**: Salinity [PSS]
- **TDS**: [mg/l]
"""
function importCOM3(path::String)
    COM3 = Dict("utc_dt" => DateTime[],
                "unix_dt" => Float64[],
                "Temp3490" => Float64[],
                "CO" => Float64[],
                "OB" => Float64[],
                "RefFuel" => Float64[],
                "TRYP" => Float64[],
                "Turb3490" => Float64[],
                "Salinity3490" => Float64[],
                "TDS" => Float64[]
                )

    for line ∈ readlines(path)
        splitline = split(line, ",")

        if length(splitline) == 13
            date = splitline[1]
            time = splitline[2]
            dt = DateTime(date*" "*time, "m/d/y H:M:S") + Year(2000)
            unix_dt = datetime2unix(dt)

            push!(COM3["utc_dt"], dt)
            push!(COM3["unix_dt"], unix_dt)

            # Temp3488
            temp = parse(Float64, splitline[6])
            push!(COM3["Temp3490"], temp)

            # CO
            co = parse(Float64, splitline[7])
            push!(COM3["CO"], co)

            # OB
            ob = parse(Float64, splitline[8])
            push!(COM3["OB"], ob)

            # RefFuel
            refuel = parse(Float64, splitline[9])
            push!(COM3["RefFuel"], refuel)

            # TRYP
            tryp = parse(Float64, splitline[10])
            push!(COM3["TRYP"], tryp)

            #Turb3490
            turb = parse(Float64, splitline[11])
            push!(COM3["Turb3490"], turb)

            #Salinity3490
            sal = parse(Float64, splitline[12])
            push!(COM3["Salinity3490"], sal)

            #TDS
            tds = parse(Float64, splitline[13])
            push!(COM3["TDS"], tds)

        end
    end
    return DataFrame(COM3)
end



"""
    importLISST(path::String)

Read an LISST file into a DataFrame

- **utc_dt**: Datetime
- **unix_dt**: unix epoch time [seconds]
- **SSC**: []
"""
function importLISST(path::String)
    LISST = Dict("utc_dt" => DateTime[],
                 "unix_dt" => Float64[],
                "SSC" => Float64[]
                )

    for line ∈ readlines(path)
        splitline = split(line, ",")
        if length(splitline) == 3 && !occursin("\0", line)
            date = splitline[1]
            time = splitline[2]
            dt = DateTime(date*" "*time, "m/d/y H:M:S") + Year(2000)
            unix_dt = datetime2unix(dt)

            push!(LISST["utc_dt"], dt)
            push!(LISST["unix_dt"], unix_dt)

            # SSC
            ssc = parse(Float64, splitline[3])
            push!(LISST["SSC"], ssc)
        end
    end
    return DataFrame(LISST)
end





"""
    processBoatFiles(basepath::String, outpath::String)

Parse NMEA files into CSVs.
"""
function processBoatFiles(basepath::String, outpath::String)
    for (root, dirs, files) in walkdir(basepath)
        @showprogress for file in files
            if !(occursin("fixed", file))
                if occursin("AirMar", file)
                    name = split(file, "_")[2]
                    airmar_gps, airmar_speed = importAirMar(joinpath(root, file))
                    CSV.write(joinpath(outpath, name*"_airmar_gps.csv"), airmar_gps)
                    CSV.write(joinpath(outpath, name*"_airmar_speed.csv"), airmar_speed)
                elseif occursin("COM1", file)
                    name = split(file, "_")[2]
                    COM1 = importCOM1(joinpath(root, file))
                    CSV.write(joinpath(outpath, name*"_COM1.csv"), COM1)

                elseif occursin("COM2", file)
                    name = split(file, "_")[2]
                    COM2 = importCOM2(joinpath(root, file))
                    CSV.write(joinpath(outpath, name*"_COM2.csv"), COM2)

                elseif occursin("COM3", file)
                    name = split(file, "_")[2]
                    COM3 = importCOM3(joinpath(root, file))
                    CSV.write(joinpath(outpath, name*"_COM3.csv"), COM3)

                elseif occursin("LISST", file)
                    name = split(file, "_")[2]
                    LISST = importLISST(joinpath(root, file))
                    CSV.write(joinpath(outpath, name*"_LISST.csv"), LISST)

                elseif occursin("nmea", file) || occursin("NMEA", file)
                    name = split(file, "_")[2]
                    nmea = importNMEA(joinpath(root, file))
                    CSV.write(joinpath(outpath, name*"_nmea.csv"), nmea)
                elseif occursin("nmea", file)
                    # skip for now
                    continue
                end
            end
        end
    end
end





function combine_boat_dfs(dirpath,
                          outpath;
                          w = -97.717472,
                          n = 33.703572,
                          s = 33.700797,
                          e = -97.712413,
                          Δx=0.10
                          )



    airmar_gps_dfs = []
    airmar_speed_dfs = []
    com1_dfs = []
    com2_dfs = []
    com3_dfs = []
    lisst_dfs = []
    nmea_dfs = []

    @showprogress for f ∈ readdir(dirpath)
        if endswith(f, "gps.csv")
            push!(airmar_gps_dfs, DataFrame(CSV.File(joinpath(dirpath, f))))
        elseif endswith(f, "speed.csv")
            push!(airmar_speed_dfs, DataFrame(CSV.File(joinpath(dirpath, f))))
        elseif endswith(f, "COM1.csv")
            push!(com1_dfs, DataFrame(CSV.File(joinpath(dirpath, f))))
        elseif endswith(f, "COM2.csv")
            push!(com2_dfs, DataFrame(CSV.File(joinpath(dirpath, f))))
        elseif endswith(f, "COM3.csv")
            push!(com3_dfs, DataFrame(CSV.File(joinpath(dirpath, f))))
        elseif endswith(f, "LISST.csv")
            push!(lisst_dfs, DataFrame(CSV.File(joinpath(dirpath, f))))
        elseif endswith(f, "nmea.csv")
            push!(nmea_dfs, DataFrame(CSV.File(joinpath(dirpath, f))))
        end
    end

    airmar_gps_df = vcat(airmar_gps_dfs...)
    airmar_speed_df = vcat(airmar_speed_dfs...)
    com1_df = vcat(com1_dfs...)
    com2_df = vcat(com2_dfs...)
    com3_df = vcat(com3_dfs...)
    lisst_df = vcat(lisst_dfs...)
    nmea_df = vcat(nmea_dfs...)


    sort!(airmar_gps_df, :utc_dt)
    unique!(airmar_gps_df)

    sort!(airmar_speed_df, :utc_dt)
    unique!(airmar_speed_df)

    sort!(com1_df, :utc_dt)
    unique!(com1_df)

    sort!(com2_df, :utc_dt)
    unique!(com2_df)

    sort!(com3_df, :utc_dt)
    unique!(com3_df)

    sort!(lisst_df, :utc_dt)
    unique!(lisst_df)

    sort!(nmea_df, :utc_dt)
    nmea_df = unique(nmea_df)
    # nmea sensor has duplicate readings sometimes
    nmea_df = combine(first, groupby(sort(nmea_df, :utc_dt), :utc_dt))

    # make sure we're in the correct bounding box
    nmea_df = nmea_df[(nmea_df.latitude .> s) .& (nmea_df.latitude .< n) .& (nmea_df.longitude .> w) .& (nmea_df.longitude .< e), :]

    # filter to times within boat GPS values
    tstart = nmea_df.utc_dt[1]
    tend = nmea_df.utc_dt[end]

    println("\tcom1 df size: ", size(com1_df))
    println("\tcom2 df size: ", size(com2_df))
    println("\tcom3 df size: ", size(com3_df))
    println("\tairmar df size: ", size(airmar_speed_df))
    println("\tlisst df size: ", size(lisst_df))


    # filter df's to those times that fall within Boat GPS times (i.e. nmea)
    com1_df_filtered = com1_df[(com1_df.utc_dt .>= tstart) .& (com1_df.utc_dt .<= tend ), :];
    com2_df_filtered = com2_df[(com2_df.utc_dt .>= tstart) .& (com2_df.utc_dt .<= tend ), :];
    com3_df_filtered = com3_df[(com3_df.utc_dt .>= tstart) .& (com3_df.utc_dt .<= tend ), :];
    airmar_speed_df_filtered = airmar_speed_df[(airmar_speed_df.utc_dt .>= tstart) .& (airmar_speed_df.utc_dt .<= tend ), :];
    lisst_df_filtered = lisst_df[(lisst_df.utc_dt .>= tstart) .& (lisst_df.utc_dt .<= tend ), :];

    println("\tcom1 filtered df size: ", size(com1_df_filtered))
    println("\tcom2 filtered df size: ", size(com2_df_filtered))
    println("\tcom3 filtered df size: ", size(com3_df_filtered))
    println("\tairmar filtered df size: ", size(airmar_speed_df_filtered))
    println("\tlisst filtered df size: ", size(lisst_df_filtered))


    # now let's interpolate to our GPS times
    interpolated = Dict()

    interpolated["longitude"] = nmea_df.longitude
    interpolated["latitude"] = nmea_df.latitude

    interpolated["unix_dt"] = nmea_df.unix_dt
    interpolated["utc_dt"] = nmea_df.utc_dt



    # go through COM1
    println("\tInterpolating COM1")
    names(com1_df_filtered)
    Br_interp = CubicSpline(com1_df_filtered.Br, com1_df_filtered.unix_dt)
    interpolated["Br"] = Br_interp.(interpolated["unix_dt"])

    Ca_interp = CubicSpline(com1_df_filtered.Ca, com1_df_filtered.unix_dt)
    interpolated["Ca"] = Ca_interp.(interpolated["unix_dt"])

    Cl_interp = CubicSpline(com1_df_filtered.Cl, com1_df_filtered.unix_dt)
    interpolated["Cl"] = Cl_interp.(interpolated["unix_dt"])

    HDO_interp = CubicSpline(com1_df_filtered.HDO, com1_df_filtered.unix_dt)
    interpolated["HDO"] = HDO_interp.(interpolated["unix_dt"])

    HDO_percent_interp = CubicSpline(com1_df_filtered.HDO_percent, com1_df_filtered.unix_dt)
    interpolated["HDO_percent"] = HDO_percent_interp.(interpolated["unix_dt"])

    NH4_interp = CubicSpline(com1_df_filtered.NH4, com1_df_filtered.unix_dt)
    interpolated["NH4"] = NH4_interp.(interpolated["unix_dt"])

    NO3_interp = CubicSpline(com1_df_filtered.NO3, com1_df_filtered.unix_dt)
    interpolated["NO3"] = NO3_interp.(interpolated["unix_dt"])

    Na_interp = CubicSpline(com1_df_filtered.Na, com1_df_filtered.unix_dt)
    interpolated["Na"] = Na_interp.(interpolated["unix_dt"])

    Salinity3488_interp = CubicSpline(com1_df_filtered.Salinity3488, com1_df_filtered.unix_dt)
    interpolated["Salinity3488"] = Salinity3488_interp.(interpolated["unix_dt"])

    SpCond_interp = CubicSpline(com1_df_filtered.SpCond, com1_df_filtered.unix_dt)
    interpolated["SpCond"] = SpCond_interp.(interpolated["unix_dt"])

    TDS_interp = CubicSpline(com1_df_filtered.TDS, com1_df_filtered.unix_dt)
    interpolated["TDS"] = TDS_interp.(interpolated["unix_dt"])

    Temp3488_interp = CubicSpline(com1_df_filtered.Temp3488, com1_df_filtered.unix_dt)
    interpolated["Temp3488"] = Temp3488_interp.(interpolated["unix_dt"])

    Turb3488_interp = CubicSpline(com1_df_filtered.Turb3488, com1_df_filtered.unix_dt)
    interpolated["Turb3488"] = Turb3488_interp.(interpolated["unix_dt"])

    pH_interp = CubicSpline(com1_df_filtered.pH, com1_df_filtered.unix_dt)
    interpolated["pH"] = pH_interp.(interpolated["unix_dt"])

    pH_mV_interp = CubicSpline(com1_df_filtered.pH_mV, com1_df_filtered.unix_dt)
    interpolated["pH_mV"] = pH_mV_interp.(interpolated["unix_dt"])


    # go through COM2
    println("\tInterpolating COM2")
    names(com2_df_filtered)
    CDOM_interp = CubicSpline(com2_df_filtered.CDOM, com2_df_filtered.unix_dt)
    interpolated["CDOM"] = CDOM_interp.(interpolated["unix_dt"])

    Chl_interp = CubicSpline(com2_df_filtered.Chl, com2_df_filtered.unix_dt)
    interpolated["Chl"] = Chl_interp.(interpolated["unix_dt"])

    ChlRed_interp = CubicSpline(com2_df_filtered.ChlRed, com2_df_filtered.unix_dt)
    interpolated["ChlRed"] = ChlRed_interp.(interpolated["unix_dt"])

    Temp3489_interp = CubicSpline(com2_df_filtered.Temp3489, com2_df_filtered.unix_dt)
    interpolated["Temp3489"] = Temp3489_interp.(interpolated["unix_dt"])

    Turb3489_interp = CubicSpline(com2_df_filtered.Turb3489, com2_df_filtered.unix_dt)
    interpolated["Turb3489"] = Turb3489_interp.(interpolated["unix_dt"])

    bg_interp = CubicSpline(com2_df_filtered.bg, com2_df_filtered.unix_dt)
    interpolated["bg"] = bg_interp.(interpolated["unix_dt"])

    bgm_interp = CubicSpline(com2_df_filtered.bgm, com2_df_filtered.unix_dt)
    interpolated["bgm"] = bgm_interp.(interpolated["unix_dt"])


    # go through COM3
    println("\tInterpolating COM3")
    names(com3_df_filtered)
    CO_interp = CubicSpline(com3_df_filtered.CO, com3_df_filtered.unix_dt)
    interpolated["CO"] = CO_interp.(interpolated["unix_dt"])

    OB_interp = CubicSpline(com3_df_filtered.OB, com3_df_filtered.unix_dt)
    interpolated["OB"] = OB_interp.(interpolated["unix_dt"])

    RefFuel_interp = CubicSpline(com3_df_filtered.RefFuel, com3_df_filtered.unix_dt)
    interpolated["RefFuel"] = RefFuel_interp.(interpolated["unix_dt"])

    Salinity3490_interp = CubicSpline(com3_df_filtered.Salinity3490, com3_df_filtered.unix_dt)
    interpolated["Salinity3490"] = Salinity3490_interp.(interpolated["unix_dt"])

    TDS_interp = CubicSpline(com3_df_filtered.TDS, com3_df_filtered.unix_dt)
    interpolated["TDS"] = TDS_interp.(interpolated["unix_dt"])

    TRYP_interp = CubicSpline(com3_df_filtered.TRYP, com3_df_filtered.unix_dt)
    interpolated["TRYP"] = TRYP_interp.(interpolated["unix_dt"])

    Temp3490_interp = CubicSpline(com3_df_filtered.Temp3490, com3_df_filtered.unix_dt)
    interpolated["Temp3490"] = Temp3490_interp.(interpolated["unix_dt"])

    Turb3490_interp = CubicSpline(com3_df_filtered.Turb3490, com3_df_filtered.unix_dt)
    interpolated["Turb3490"] = Turb3490_interp.(interpolated["unix_dt"])


    # go through lisst
    println("\tInterpolating LISST")
    names(lisst_df_filtered)
    # cubic spline failing for some reason.
    SSC_interp = QuadraticInterpolation(lisst_df_filtered.SSC, lisst_df_filtered.unix_dt)
    interpolated["SSC"] = []
    for t ∈ interpolated["unix_dt"]
        try
            push!(interpolated["SSC"], SSC_interp(t))
        catch e
            push!(interpolated["SSC"], NaN)
        end
    end


    # generate UTMz coords rounded to nearest Δx
    xs = zeros(length(interpolated["longitude"]))
    ys = zeros(length(interpolated["longitude"]))
    zones = []
    isnorths = []

    for i ∈ 1:length(xs)
        utmz = UTMZ(LLA(interpolated["latitude"][i], interpolated["longitude"][i]), wgs84)
        xs[i] = bump_to_nearest_Δx(utmz.x, Δx)
        ys[i] = bump_to_nearest_Δx(utmz.y, Δx)
        push!(zones, utmz.zone)
        push!(isnorths, utmz.isnorth)

        # update the lat/lon to reflect new resolution
        lla = LLAfromUTMZ(wgs84)(UTMZ(xs[i], ys[i], 0.0, zones[i], isnorths[i]))
        interpolated["latitude"][i] = lla.lat
        interpolated["longitude"][i] = lla.lon
    end

    interpolated["X"] = bump_to_nearest_Δx.(xs, Δx)
    interpolated["Y"] = bump_to_nearest_Δx.(ys, Δx)
    interpolated["zone"] = zones
    interpolated["isnorth"] = isnorths

    CSV.write(joinpath(outpath, "Targets.csv"), DataFrame(interpolated))
end






end
