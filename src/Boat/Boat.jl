module Boat


using CSV, DataFrames
using Dates

# we will need to parse NMEA encoded gps/sensor strings
# see this site: http://aprs.gids.nl/nmea/


export importAirMar, importCOM1, importCOM2, importCOM3, importLISST, importNMEA




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






end
