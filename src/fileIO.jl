export setup_examples

function setup_examples(workingdir,filetype="all")
    if filetype == "all" || filetype == "notebook" 
        if !isdir(workingdir,"examples") 
            mkdir(joinpath(workingdir,"examples"))
        end

        exdir = joinpath(dirname(pathof(ComputationalHeatTransfer)),"../examples")
        for (root, dirs, files) in walkdir(exdir)
            for file in files
                cp(joinpath(root, file),joinpath(workingdir,"examples",file),force=true)
                chmod(joinpath(workingdir,"examples",file),0o644)
            end
        end
    end

    if filetype == "all" || filetype == "experiment" 
        if !isdir(workingdir,"expdata") 
            mkdir(joinpath(workingdir,"expdata"))
        end

        exdir = joinpath(dirname(pathof(ComputationalHeatTransfer)),"../expdata")
        for (root, dirs, files) in walkdir(exdir)
            for file in files
                cp(joinpath(root, file),joinpath(workingdir,"expdata",file),force=true)
                chmod(joinpath(workingdir,"expdata",file),0o644)
            end
        end
    end    

    if filetype == "all" || filetype == "simulation" 
        if !isdir(workingdir,"numedata") 
            mkdir(joinpath(workingdir,"numedata"))
        end
    end    

    return
end
