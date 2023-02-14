export inch,gravity

const inch = 2.54e-2; 
const gravity = 9.8;


include("intube/Systems.jl")
include("intube/Thermomodel.jl")
include("intube/Tools.jl")
include("intube/Mapping.jl")
include("intube/Timemarching.jl")
include("intube/Postprocessing.jl")
include("intube/CoolProp.jl")
include("intube/Plotrecipe.jl")
include("intube/DrawingOHP.jl")
include("intube/Preprocessing.jl")
# include("intube/FluidProperty.jl")
include("intube/Heatercondenser.jl")
include("intube/callback/boiling.jl")
include("intube/callback/vapormerging.jl")
include("intube/callback/liquidmerging.jl")
include("intube/callback/fixdx.jl")
