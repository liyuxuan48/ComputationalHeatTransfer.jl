# ASETS-II cases simulation

This is an example of a simulation package for conjugate heat transfer of an oscillating heat pipe. SI units are used and units are emitted

## What do we need to solve an OHP problem?

**specify properties** : Solid property, Fluid property

**set the geometries** : Computational domain, Heaters/Condensers, OHP shapes

**construct the systems** : Fluid system(1D), HeatConduction system(2D)

**initialize** : initialize the integrators and the data structs for saving

**solve**

**save/examine**

## Packages

Firstly, let's import the necessary packages, you may need to install them for the first time.


```julia
using ComputationalHeatTransfer # our main package
using ProgressMeter # to have a progessbar when runing the simulation
using Plots # for plotting
gr()  #ploting backend (the fastest one)
```

# Specify properties

### Solid Physical parameters

params is the HeatConductionParameters for the plate material. The numbers below represents aluminum.


```julia
ρₛ = 2730; # density
cₛ  = 8.93e02; # specific heat
kₛ  = 1.93e02; # heat conductivity
plate_d = 1.5e-3; # effective d (The thickness of an ideal uniform thickness plate occupying the same volume)
params = HeatConductionParameters(ρₛ ,cₛ ,kₛ ,thickness=plate_d)
```

### Fluid Physical parameters

p_fluid contains the vapor and liquid properties at a constant reference temperature. Noted that the vapor pressure and the vapor density will be functions of temperatures during the simulation, other properties are extracted from p_fluid as an approximate value.


```julia
Tref = 291.2 # reference temperature
fluid_type = "Butane"
p_fluid = SaturationFluidProperty(fluid_type,Tref)
```

# Set the geometries

### Geometry parameters
The 2D domain is of rectangular shape (slightly different from ASETS-II). In the future it can be of arbitrary shape using the immersedlayers.jl package.


```julia
Lx = 0.1524; # plate size x [m]
Ly = 0.0648; # plate size y [m]
xlim = (-Lx/2,Lx/2) # plate x limits
ylim = (-Ly/2,Ly/2) # plate y limits
```

### Set mesh size and maximum time step for plate heat conduction
Δx is controlled by Δx = α*gridPe and set having the same order of magitute of tube diameter 1e-3. Fourier number is used to give a safety "cap" of time step you can choose in the fluid module


```julia
Δx,Δt_max = setstepsizes(params.α,gridPe=8.0,fourier=0.3)
```

### Set up the evaporators and condensers
Right now, the OHPtype looks up a preset dictionary of OHP evaporators and condensers.

You can also customize them in the OHP DIY notebook




```julia
OHPtype = "ASETS-II OHP 2 LARGE HEATER"
power = 40 # total heater power in watts
Tc = Tref; # condenser temperature
eparams,cparams = OHPConfiguration(OHPtype,power,Tc,Δx);
```

### Set up OHP channel's shape
construct_ohp_curve is a built-in function that generates two arrays: x that contains all x values of the discrete points, and y contains all y values. x and y have the same length. 

You can also customize this function to generate an OHP shape of your choice as long as they produce x array and y array.


```julia
x, y = construct_ohp_curve("ASETS",Δx) # get x and y coordinates for the channel
ohp = BasicBody(x,y) # build a BasicBody based on x,y

ohpgeom = ComputationalHeatTransfer.LineSourceParams(ohp) # build a line heat source based on BasicBody
```

### Plot what you got so far
This is a exmaple of the compuational domain (the box) and the OHP channel serpentine (in blue)


```julia
# plot ohp
plt = plot(ohp,fillalpha=0,linecolor=:black,xlims=xlim,ylims=ylim,framestyle = :box,xlabel="x [m]",ylabel="y [m]")

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

# Construct the systems

### Create HeatConduction system
The solid module dealing with the 2D conduction, evaporator, condenser, and the OHP line heat source is constructed here.


```julia
sys_plate = HeatConduction(params,Δx,xlim,ylim,Δt_max,qline=ohpgeom,qflux=eparams,qmodel=cparams)
```

### Create OHP inner channel system
sys_tube: fluid module system


```julia
sys_tube = initialize_ohpsys(OHPtype,fluid_type,sys_plate,p_fluid,Tref,power)
```

# Initialize

### set time step


```julia
tspan = (0.0, 2.0); # start time and end time
dt_record = 0.01   # saving time interval

tstep = 1e-3     # actrual time marching step
```

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

# Solve

### Run the simulation and store data


```julia
@showprogress for t in tspan[1]:tstep:tspan[2]

    timemarching!(integrator_tube,integrator_plate,tstep)

    if (mod(integrator_plate.t,dt_record) < 1e-6) || (mod(-integrator_plate.t,dt_record) < 1e-6)
        store!(sr,integrator_tube,integrator_plate)
    end

end
```

# Store data


```julia
save_path = dirname(dirname(dirname(pwd())))*"/OHPnume/OHP2_40W_large.jld2"
save(save_path,"SimulationResult",sr)
```

---

*This notebook was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*
