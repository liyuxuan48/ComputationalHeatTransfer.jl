```julia
using Pkg
Pkg.activate(dirname(pwd())) # using current environment for development
```


```julia
using Interact
using UnPack
```


```julia
using ComputationalHeatTransfer
using LaTeXStrings
using JLD2
using Interpolations
using Plots
gr()  
```

# Read simulation data


```julia
read_path = dirname(dirname(dirname(pwd())))*"/OHPnume/OHP2_40W_large.jld2"
SimuResult = load(read_path)["SimulationResult"];
```

### get time array


```julia
t = SimuResult.tube_hist_t;
```

# Plot 2D graphs

### film and slug dynamics


```julia
@gif for i in eachindex(t)
    plot(OHPSlug(),i,SimuResult)
end
```

### plate T [K]


```julia
@gif for i in eachindex(t)
    plot(OHPTemp(),i,SimuResult,clim=(291.2,292.0))
end
```

### 2D superheat


```julia
@gif for i in eachindex(t)
# @gif for i =1:10
    plot(OHPSuper(),i,SimuResult)
end
```

### 2D pressure


```julia
@gif for i in eachindex(t)
    plot(OHPPres(),i,SimuResult)
end
```

# Plot 2D interpolated curves

### Interpolate 2D T data from the plate for fixed sensors on the plate

place the 2D sensors


```julia
x2Dsensors = [-2.75,-1.4,-0.8,0.0,0.0,0.8,1.4,2.75] .* inches
y2Dsensors = [0.0,   0.0, 0.0,0.0,0.4,0.0,0.0,0.0] .* inches
plate_sensors = (x2Dsensors,y2Dsensors);
```

get the curve


```julia
t_hist,g_hist = getTcurve(plate_sensors,SimuResult);
```

### Read experiment T data


```julia
import XLSX
```

read experiment file (customizable)


```julia
# customize 
expfile = expfileDict["O002_H001_P040"]
exppath = dirname(dirname(dirname(pwd())))*"/OHPexp/"
xf = XLSX.readxlsx(exppath*expfile);
```

get experiment data


```julia
Onum, Hnum, power_exp = getconfig(expfile)
RTDt,RTD = getRTD(xf,Onum);
```

### 2D interpolated temperature curve at fixed sensors


```julia
RTD_for_plotting = [1,4,8];
```


```julia
# plot OHP
plot(OHP(),SimuResult)

# plot sensors
scatter!(x2Dsensors[RTD_for_plotting],y2Dsensors[RTD_for_plotting])
annotate!(x2Dsensors[RTD_for_plotting], y2Dsensors[RTD_for_plotting].-0.005, RTD_for_plotting)
```


```julia
plot(OHPTcurve(),RTD_for_plotting,(t_hist,g_hist),SimuResult)
plot!(OHPTexp() ,RTD_for_plotting,(RTDt,RTD)     ,SimuResult)
```

### 2D interpolated thermal conductance


```julia
ihot = 4 # hot sensor  for calculating thermal conductance
icold = 8 # cold sensor  for calculating thermal conductance;
```


```julia
# plot them separately
plot(OHPCond(),(ihot,icold),(t_hist,g_hist),(RTDt,RTD),SimuResult)
```

### Liquid slug velocity statistics


```julia
# fix title and ylabel and legend
plot(OHPV(), SimuResult::SimulationResult,ylimit=(-2,2))  
```

# Plot 1D interpolated curves


```julia
# tell the time
@manipulate for i in 1:1:length(t)
    plot(OHP1DT(),i,SimuResult,xlim=(1,2))
#     plot!(twinx(),OHPTwall(),i,SimuResult,xlim=(1,2))
    plot!(twinx(),OHP1DΔT(),i,SimuResult,xlim=(1,2))
#     plot!(twinx(),OHP1DP(),i,SimuResult,xlim=(1,2))
end
```

### 1D sensor selector


```julia
L = SimuResult.integrator_tube.p.tube.L
@manipulate for ξ in 0:1e-3:L
    plot(OHP(),SimuResult) # plot the ohp layout

    xprobe,yprobe = oneDtwoDtransform(ξ,SimuResult)
    scatter!([xprobe],[yprobe])
end
```

## Plot 1D property curve for a fixed location sensor


```julia
xsensors1D = [2.097, 3.0, 4,4.1]
```


```julia
θhist1D,phist1D = get1DTandP(xsensors1D, SimuResult);
```


```julia
plot(t,θhist1D,label=string.("ξ=", xsensors1D'),xlabel="time [s]", ylabel="temperature [K]")
```

### get boiling data


```julia
boil_data,boil_num_x,boil_num_t,t_boil,x2D_boil,y2D_boil,boil_dt = get_boil_matrix(SimuResult::SimulationResult);
```


```julia
plot(OHP(),SimuResult)
scatter!(x2D_boil,y2D_boil,
    colorbar=true,markeralpha=delta.(boil_num_x),colorbar_title="\n boiling frequency [Hz]",right_margin=3Plots.mm,marker_z=boil_num_x./SimuResult.tube_hist_t[end],markerstrokewidth=0,markercolor=cgrad(:greys, rev = true))
```


```julia
plot(t_boil,boil_num_t./boil_dt,
color=:orange, legend=:topleft, ylabel="f [HZ]",xlabel="time [s]", label="overall boiling frequency")
```
