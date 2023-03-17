# DIY

This notebooks shows how to custimize the heater/condenser and ohp configuration

## Packages
If you need to update the ComputationalHeatTransfer package, I suggest you activate the package's path as your current environment. If you don't need to modify the code within the packge, you don't have to run this block.


```julia
using Pkg
Pkg.activate(dirname(pwd())) # using current environment for development
```

    [32m[1m  Activating[22m[39m project at `~/Documents/GitHub/ComputationalHeatTransfer.jl`


Firstly, let's import the necessary packages, you may need to install them for the first time.


```julia
using ComputationalHeatTransfer # our main package
using LaTeXStrings # for latex strings
using JLD2 # for file input/output
using ProgressMeter # to have a progessbar when runing the simulation
using XLSX # for reading experimental data in CSV format
using Plots # for plotting
gr()  #ploting backend (the fastest one)
```




    Plots.GRBackend()



# Properies

## Solid Physical parameters

params is the HeatConductionParameters for the plate material.

The numbers below represents aluminum.


```julia
ρₛ = 2730; # density
cₛ  = 8.93e02; # specific heat
kₛ  = 1.93e02; # heat conductivity
plate_d = 1.5e-3; # effective d (The thickness of an ideal uniform thickness plate occupying the same volume)
params = HeatConductionParameters(ρₛ ,cₛ ,kₛ ,thickness=plate_d)
```




    HeatConductionParameters(2730.0, 893.0, 193.0, 7.916682048820907e-5, 0.0015)



## Fluid Physical parameters

p_fluid contains the vapor and liquid properties at a constant reference temperature. Noted that the vapor pressure and the vapor density will be functions of temperatures during the simulation, other properties are extracted from p_fluid as an approximate value.


```julia
Tref = 291.2 # reference temperature
fluid_type = "Butane"
p_fluid = SaturationFluidProperty(fluid_type,Tref)
```




    Saturation properties for Butane at constant temperature 291.2 [K]




# Plate Conduction Part

### Geometry parameters
The 2D domain is of rectangular shape (slightly different from ASETS-II). In the future it can be of arbitrary shape using the immersedlayers.jl package.


```julia
Lx = 0.1524; # plate size x [m]
Ly = 0.0648; # plate size y [m]
xlim = (-Lx/2,Lx/2) # plate x limits
ylim = (-Ly/2,Ly/2) # plate y limits
```




    (-0.0324, 0.0324)



### Set mesh size and maximum time step for plate heat conduction
Δx is controlled by Δx = α*gridPe and set having the same order of magitute of tube diameter 1e-3. Fourier number is used to give a safety "cap" of time step you can choose in the fluid module


```julia
Δx,Δt_max = setstepsizes(params.α,gridPe=8.0,fourier=0.3)
```




    (0.0006333345639056725, 0.001520002953373614)



### Set up the evaporators and condensers
Right now, the OHPtype looks up a preset dictionary of OHP evaporators and condensers.

You can also customize them in the OHP DIY notebook.

Firstly let's give the total heater power


```julia
power = 10 # total heater power in watts
```




    10



Then let's construct a heater


```julia
Lheater_x = Lx*0.1
Lheater_y = Ly*0.9

qe = power/Lheater_x/Lheater_y

eb1 = Rectangle(Lheater_x/2,Lheater_x/2,1.5*Δx)
Tfe = RigidTransform((0.0,-0.0),0.0)
Tfe(eb1)

eparams = [PrescribedHeatFluxRegion(qe,eb1)];
```

Then let's consctruct a condenser


```julia
Lcondenser_x = Lx*0.2
Lcondenser_y = Ly*0.9

hc = 2000.0
qe = power/Lheater_x/Lheater_y

cb1 = Rectangle(Lheater_y/2,Lheater_y/2,1.5*Δx)
Tfc = RigidTransform((Lx*0.3,-0.0),0.0)
Tfc(cb1)

Tc = Tref
cparams = [PrescribedHeatModelRegion(hc,Tc,cb1)];
```

### Set up OHP channel's shape
in this example we construct a user-defined curve


```julia
a = 0.03
θ = 0:2π/1000:2π
# r = a*(1 .- sin.(θ))
r = a*sin.(4θ)
x = r .* cos.(θ)
y = r .* sin.(θ);
```


```julia
# x, y = construct_ohp_curve("ASETS",Δx) # get x and y coordinates for the channel
ohp = BasicBody(x,y) # build a BasicBody based on x,y
```




    Basic pointwise-specified body with 1001 points
       Current position: (0.0,0.0)
       Current angle (rad): 0.0




### Plot what you got so far
This is a exmaple of the compuational domain (the box) and the OHP channel serpentine (in blue)


```julia
# plot ohp
plt = plot(ohp,fillalpha=0,linecolor=:black,xlims=xlim,ylims=ylim,framestyle = :box)

# plot heaters (red)
for ep in eparams
    plot!(ep)
end

# plot condensers (blue)
for cp in cparams
    plot!(cp)
end

# show plot
plt
```




    
![svg](output_28_0.svg)
    



### Create HeatConduction system
The solid module dealing with the 2D conduction, evaporator, condenser, and the OHP line heat source is constructed here.


```julia
ohpgeom = ComputationalHeatTransfer.LineSourceParams(ohp) # build a line heat source based on BasicBody
sys_plate = HeatConduction(params,Δx,xlim,ylim,Δt_max,qline=ohpgeom,qflux=eparams,qmodel=cparams)
```




    Unbounded Heat conduction system on a grid of size 252 x 112 and 0 static immersed points




### Create OHP inner channel system
sys_tube: fluid module system


```julia
sys_tube = initialize_ohpsys(OHPtype,fluid_type,sys_plate,p_fluid,Tref,power)
```




    1001 point OHP system filled with Butane




### set time step


```julia
tspan = (0.0, 2.0); # start time and end time
dt_record = 0.01   # saving time interval

tstep = 1e-3    # actrual time marching step
```




    0.001



### combine inner tube and plate together


```julia
u_plate = newstate(sys_plate) .+ Tref # initialize plate T field to uniform Tref
integrator_plate = init(u_plate,tspan,sys_plate) # construct integrator_plate

u_tube = newstate(sys_tube) # initialize OHP tube 
integrator_tube = init(u_tube,tspan,sys_tube); # construct integrator_tube
```

### initialize arrays for saving


```julia
sr = SimulationResult(integrator_tube,integrator_plate);
```

### Run the simulation and store data


```julia
@showprogress for t in tspan[1]:tstep:tspan[2]

    timemarching!(integrator_tube,integrator_plate,tstep)

    if (mod(integrator_plate.t,dt_record) < 1e-6) || (mod(-integrator_plate.t,dt_record) < 1e-6)
        store!(sr,integrator_tube,integrator_plate)
    end

end
```

    [32mProgress: 100%|█████████████████████████████████████████| Time: 0:00:26[39m


## Store data


```julia
save_path = dirname(dirname(dirname(pwd())))*"/OHPnume/OHP2_40W_large.jld2"
save(save_path,"SimulationResult",sr)
```

---

*This notebook was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*
