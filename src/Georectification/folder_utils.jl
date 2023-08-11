"""
    get_bil_files(dir, file_id)

Get a list of all .bil files in `dir` whose filename matches `file_id`. Returned list is sorted by the capture number at end of filename string.
"""
function get_bil_files(dir::String, file_id::String)
    bils = []
    for (root, dirs, files) ∈ walkdir(dir)
        for file ∈ files
            ffull = joinpath(root, file)
            namebase = split(split(ffull, "/")[end-1], "-")[1]

            if endswith(file, ".bil") && file_id == namebase
                push!(bils, joinpath(root, file))
            end
        end
    end

    # get list of file numbers to produce sorting indices

    endings = [split(f, "_")[end] for f ∈ bils]
    number = [lpad(split(f, "-")[1], 2, "0") for f ∈ endings]

    idx = sortperm(number)

    return bils[idx]
end



"""
    get_raw_file_list(bilpath)

Given a bil path, return the full list of paths to files needed for georectification.
"""
function get_raw_file_list(bilpath::String)
    basepath = "/"*joinpath(split(bilpath, "/")[1:end-1]...)
    bilhdrpath = ""
    timespath = ""
    specpath = ""
    spechdrpath = ""
    lcfpath = ""

    for f ∈ readdir(basepath)
        if endswith(f, ".bil.hdr")
            bilhdrpath = joinpath(basepath, f)
        elseif endswith(f, ".lcf")
            lcfpath = joinpath(basepath, f)
        elseif endswith(f, ".times")
            timespath = joinpath(basepath, f)
        elseif endswith(f, ".spec")
            specpath = joinpath(basepath, f)
        elseif endswith(f, ".spec.hdr")
            spechdrpath = joinpath(basepath, f)
        else
            continue
        end
    end

    return (;bilpath=bilpath, bilhdr=bilhdrpath, lcfpath=lcfpath, timespath=timespath, specpath=specpath, spechdr=spechdrpath)
end


