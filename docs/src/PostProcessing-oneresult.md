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

## Get fluid properties


```julia
# fluid_type = "butane"
Tᵥ = 291.2
```

## Read data


```julia
read_path = dirname(dirname(dirname(pwd())))*"/OHPnume/OHP2_40W_large.jld2"
# read_path = dirname(dirname(dirname(pwd())))*"/Hoffman/ohp82/OHPnume/plus_0.2_minus_0.19_hc2ratio_45.jld2"
OHPdata = load(read_path);
```


```julia
sysfinal,ohp,integrator_plate,integrator_tube,boil_hist,plate_T_hist = translateOHPdata(OHPdata);
```


```julia
x = [-2.75,-1.4,-0.8,0.0,0.0,0.8,1.4,2.75] .* inches
y = [0.0,   0.0, 0.0,0.0,0.4,0.0,0.0,0.0] .* inches
ghist,thist = getTcurve(x,y,OHPdata);
```

## Get experiment RTD result


```julia
import XLSX
```


```julia
expfile = expfileDict["O002_H001_P040"]
exppath = dirname(dirname(dirname(pwd())))*"/OHPexp/"
xf = XLSX.readxlsx(exppath*expfile)
Onum, Hnum, power_exp = getconfig(expfile)
RTD,RTDt = getRTD(xf,Onum);
```

## boiling data


```julia
boil_data,boil_num_x,boil_num_t,t_boil,x2D_boil,y2D_boil = get_boil_matrix(OHPdata);
# scatter(boil_matrix)
```

## 2D graphs

### film and dynamics


```julia
gr()
Hₗ = sysfinal[1].liquid.Hₗ
adjust = 1e-2;

anim = @animate for i=1:length(sysfinal)
    Htmp = sys_to_Harray(sysfinal[i])
    Htmp_marker = round.(div.(Htmp,Hₗ-1e-10))
    plot(ohp,clim=(0,2),fillalpha=0,linewidth=2.0,linecolor=palette([:red,  :blue, :yellow]),
        line_z=Htmp_marker,xlabel="x ",ylabel="y ",border=:none,axis=nothing,bbox_inches="tight")
        annotate!(0.0, 0.028, string("time = ", round(thist[i], digits=2), "[s]"), :black,legend=false)
    scatter!([-0.066+adjust],[-0.028],color=:yellow);scatter!([-0.03+adjust],[-0.028],color=:red);scatter!([0.02+adjust],[-0.028],color=:blue);
    annotate!(-0.05+0.002+adjust, -0.028, "vapor with film", :black)
    annotate!(-0.01+0.005+adjust, -0.028, "dry vapor", :black)
    annotate!(0.03+0.002+adjust, -0.028, "liquid", :black)
end
gif(anim, "slug_fps15.gif", fps = 30)
```

### plate T [K]


```julia
gr()
Tmax = maximum(plate_T_hist[end])
Tmin = minimum(plate_T_hist[1])
grid = OHPdata["integrator_plate"].p.grid
xlim = grid.xlim[1]
ylim = grid.xlim[2]

anim = @animate for i = 1:1:length(sysfinal)
# @gif for i = 1:1:1
heatmap(plate_T_hist[i],grid,legend=true,color=cgrad(:thermal),
        xlimit=xlim,ylimit=ylim,clim=(Tmin,Tmax),line_z=0,xlabel="x [m]",ylabel="y [m]",
        colorbar_title = "\n T[K]",right_margin = 5Plots.mm)
scatter!([x[1],x[4],x[8]],[y[1],y[4],y[8]])
annotate!(x[1]+0.002, y[1]+0.005, "RTD1", :white)
annotate!(x[4]+0.002, y[4]+0.005, "RTD4", :white)
annotate!(x[8]-0.001, y[8]+0.005, "RTD8", :white)
annotate!(0.05, -0.028, string("time = ", round(thist[i], digits=2), "[s]"), :white,legend=false)
end
gif(anim, "temperature_fps15.gif", fps = 30)
```

## 1D curve

## temperature curve


```julia
Tᵥ = 291.2
```


```julia
label_for_plotting = [1,4,8];
```


```julia
p1 = plot()
plot_title = "temperature curve: OHP2 large heater"
for index in eachindex(label_for_plotting)
    i = label_for_plotting[index]
    if index == 1
        p1 = plot(thist,ghist[i] .-Tᵥ, color=palette(:tab10)[index], label=string("RTD", i," simulation"),linewidth=2,legend = :topleft)
        scatter!(RTDt .- RTDt[1],RTD[:,i] .- RTD[1,i],color=palette(:tab10)[index],label=string("RTD", i," experiment"))
    else
        plot!(thist,ghist[i] .-Tᵥ,color=palette(:tab10)[index],label=string("RTD", i," simulation"),linewidth=2)
        scatter!(RTDt .- RTDt[1],RTD[:,i] .- RTD[1,i],color=palette(:tab10)[index], label=string("RTD", i," experiment"),
            xlabel="time [s]",ylabel="T-T₀ [K]",xlim=(0,thist[end]),ylim=(0,60),title=(plot_title))
    end
end
p1
```


```julia
gr()
boil_dt = 0.1
i1 = 4 #RTD number
i2 = 8 #RTD number

power = 40

p1 = plot(thist,power./(ghist[i1] .-ghist[i2]), right_margin=10Plots.mm,label=string("conductance simulation"),linewidth=2,legend=:topright)
scatter!(RTDt .- RTDt[1],power./(RTD[:,i1] .- RTD[:,i2]), label=string("conductance experiment"),
    xlim=(0,300),ylim=(0,15),title="thermal conductance and boiling frequency",xlabel="time [s]",ylabel="P/(T₄-T₈) [W/K]")
plot!(twinx(), t_boil,boil_num_t./boil_dt, color=:orange, legend=:topleft, ylabel="f [HZ]",ylim=(-2,300),xlim=(0,thist[end]),label="overall boiling frequency")
plot!(right_bottom=10Plots.mm)
# savefig(p1,"cond_curve normal.pdf")
```

## dXdt distribution


```julia
using Statistics
using EasyFit
```


```julia
Vavg_hist = []
Vmax_hist = []
Vmin_hist = []
for sysi in sysfinal
    V = [elem[2] for elem in sysi.liquid.dXdt]
#     Vavg = mean(abs.(V))
    Vavg = mean(V)
    Vmax = maximum((V))
    Vmin = minimum((V))
    
    push!(Vavg_hist, Vavg)
    push!(Vmax_hist, Vmax)
    push!(Vmin_hist, Vmin)
end
```


```julia
plot(thist,movavg(Float64.(Vmax_hist),10).x,label="v range",fillalpha = 0.35, c = :blue,fillrange = movavg(Float64.(Vmin_hist),10).x,fillcolor=:green)
plot!(thist,movavg(Float64.(Vavg_hist),10).x,label="v avg")
plot!(thist,movavg(Float64.(Vmin_hist),10).x,ylim=(-2,2),color=:blue, label=false, xlabel="t [s]", ylabel="v [m/s]",title="velocity distribution")
# savefig("velocity dryout.pdf")
```


```julia
U = mean(Vavg_hist)
```

## Plot 1D property snapshots


```julia
plot(sysfinal[end],plottype="P")
```

## sensor selector


```julia
grid = OHPdata["integrator_plate"].p.grid
xlim = grid.xlim[1]
ylim = grid.xlim[2]

tube_sys = OHPdata["integrator_tube"].p
L = tube_sys.tube.L
@manipulate for ξ in 0:1e-3:L
    plot(ohp,fillalpha=0,linecolor=:blue,xlims=xlim,ylims=ylim,framestyle = :box)
    x2,y2 = oneDtwoDtransform([ξ],OHPdata)
#     println(typeof(x2))
    scatter!([x2],[y2])
end
```

## Plot 1D property curve for a fixed location sensor


```julia
xsensors = [1.0,1.5]
```


```julia
θhist = []
phist = []
for systemp in sysfinal
    ptemp = systemp.mapping.P_interp_liquidtowall[xsensors]
    θtemp = PtoT.(ptemp)
    push!(phist, ptemp)
    push!(θhist, θtemp)
end
phist = hcat(phist...)'
θhist = hcat(θhist...)';
```


```julia
plot(θhist)
```


```julia
plot(phist)
```


```julia

```
