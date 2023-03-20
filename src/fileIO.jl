export setup_examples

function setup_examples(workingdir,filetype="all")
    if filetype == "all" || filetype == "notebook" 
        exdir = joinpath(dirname(pathof(ComputationalHeatTransfer)),"../examples")
        for (root, dirs, files) in walkdir(exdir)
            for file in files
                cp(joinpath(root, file),joinpath(workingdir,file),force=true)
                chmod(joinpath(workingdir,file),0o644)
            end
        end
    end

    if filetype == "all" || filetype == "experiment" 
        exdir = joinpath(dirname(pathof(ComputationalHeatTransfer)),"../expdata")
        for (root, dirs, files) in walkdir(exdir)
            for file in files
                cp(joinpath(root, file),joinpath(workingdir,file),force=true)
                chmod(joinpath(workingdir,file),0o644)
            end
        end
    end    
end