#   # DIY

#   This notebooks shows how to custimize the heater/condenser and ohp configuration

#   # Packages

#   Firstly, let's import the necessary packages, you may need to install them
#   for the first time.

using ComputationalHeatTransfer # our main package
using Plots # for plotting
using ProgressMeter # to have a progress bar in the calculation

#   # Specify properties

#   ### Solid Physical parameters

#   params is the HeatConductionParameters for the plate material. The numbers
#   below represents aluminum.

ρₛ = 2730; # material density [kg/m^3]
cₛ  = 8.93e02; # material specific heat [J/kg K]
kₛ  = 1.93e02; # material heat conductivity
plate_d = 1.5e-3; # effective d (The thickness of an ideal uniform thickness plate occupying the same volume)
params = HeatConductionParameters(ρₛ ,cₛ ,kₛ ,thickness=plate_d)

#   ### Fluid Physical parameters

#   pfluid contains the vapor and liquid properties at a constant reference
#   temperature. Noted that the vapor pressure and the vapor density will be
#   functions of temperatures during the simulation, other properties are
#   extracted from pfluid as an approximate value.

Tref = 291.2 # reference temperature
fluid_type = "Butane"
p_fluid = SaturationFluidProperty(fluid_type,Tref)

#   # Set the geometries

#   ### Geometry parameters

#   The 2D domain is of rectangular shape (slightly different from ASETS-II). In
#   the future it can be of arbitrary shape using the immersedlayers.jl package.

Lx = 0.1524; # plate size x [m]
Ly = 0.0648; # plate size y [m]
xlim = (-Lx/2,Lx/2) # plate x limits
ylim = (-Ly/2,Ly/2) # plate y limits

#   ### Set mesh size and maximum time step for plate heat conduction

#   Δx is controlled by Δx = α*gridPe and set having the same order of magitute
#   of tube diameter 1e-3. Fourier number is used to give a safety "cap" of time
#   step you can choose in the fluid module

Δx,Δt_max = setstepsizes(params.α,gridPe=8.0,fourier=0.3)

#   # Set up the evaporators and condensers

#   In the "OHP simulation" notebook, I use "OHPtype" to look up a preset dictionary of OHP evaporators and condensers.

#   You can also customize them, following the procedure below in this notebook.

#   Firstly let's give the total heater power

power = 30 # total heater power in watts

# Then let's construct a heater

Lheater_x = Lx*0.1
Lheater_y = Ly*0.9

qe = power/Lheater_x/Lheater_y

eb1 = Rectangle(Lheater_x/2,Lheater_x/2,1.5*Δx)
Tfe = RigidTransform((-Lx*0.1,Ly*0.1),0.0)
Tfe(eb1)

eparams = [PrescribedHeatFluxRegion(qe,eb1)];

# Then let's consctruct a condenser

Lcondenser_x = Lx*0.2
Lcondenser_y = Ly*0.9

hc = 2000.0 # condenser heat transfer coefficient

cb1 = Rectangle(Lheater_y/2,Lheater_y/2,1.5*Δx)
Tfc = RigidTransform((Lx*0.3,-0.0),0.0)
Tfc(cb1)

Tc = Tref
cparams = [PrescribedHeatModelRegion(hc,Tc,cb1)];

# # Set up OHP channel's shape

# Similarly, In the "OHP simulation" notebook, I used **construct_ohp_curve("ASETS",Δx)** to look up a preset dictionary of ASETS-II OHP.

# You can customize the ohp curve in either of the two ways:

#  1. simply supply two arrays of x and y of the same length:

a = 0.03
θ = 0:2π/1000:2π
r = a*sin.(2θ)
x = r .* cos.(θ)
y = r .* sin.(θ);

plot(x,y,aspectratio=1)

# 2. **construct_ohp_curve(nturn, pitch, height, gap, ds, x0, y0, flipx, flipy, angle)**, a built-in function to generate a closed loop multi-turn channel

ds = 1.5Δx # point interval
nturn = 6 # number of turns
width_ohp = 30*1e-3 
length_ohp = 70*1e-3
gap = 1e-3 # gap between the closed loop end to the channel(not the distance between each channels)
pitch = width_ohp/(2*nturn+1) # pitch between channels
rotation_angle = 3π/8
x0, y0 = -length_ohp/2 * 1.02, -width_ohp/2 * 0.1 # starting point location

x,y = construct_ohp_curve(nturn,pitch,length_ohp,gap,ds,x0,y0,false,false,rotation_angle)

plot(x,y,aspectratio=1)

# 
ohp = BasicBody(x,y) # build a BasicBody based on x,y
ohpgeom = ComputationalHeatTransfer.LineSourceParams(ohp) # build a line heat source based on BasicBody

#   ### Plot what you got so far

#   This is a exmaple of the compuational domain (the box) and the OHP channel
#   serpentine (in blue)

## plot ohp
plt = plot(ohp,fillalpha=0,linecolor=:black,xlims=xlim,ylims=ylim,framestyle = :box,xlabel="x [m]",ylabel="y [m]")

## plot heaters (red)
for ep in eparams
    plot!(ep)
end

## plot condensers (blue)
for cp in cparams
    plot!(cp)
end

## show plot
plt

#   # Construct the systems

#   ### Create HeatConduction system

#   The solid module dealing with the 2D conduction, evaporator, condenser, and
#   the OHP line heat source is constructed here.

sys_plate = HeatConduction(params,Δx,xlim,ylim,Δt_max,qline=ohpgeom,qflux=eparams,qmodel=cparams)

#   ### Create OHP inner channel system

#   sys_tube: fluid module system

sys_tube = initialize_ohpsys(sys_plate,p_fluid,power)

#   # Initialize

#   ### set time step


tspan = (0.0, 5.0); # start time and end time
dt_record = 0.01   # saving time interval

tstep = 1e-3     # actrual time marching step

#   ### combine inner tube and plate together

u_plate = newstate(sys_plate) .+ Tref # initialize plate T field to uniform Tref
integrator_plate = init(u_plate,tspan,sys_plate) # construct integrator_plate

u_tube = newstate(sys_tube) # initialize OHP tube 
integrator_tube = init(u_tube,tspan,sys_tube); # construct integrator_tube

#   ### initialize arrays for saving


SimuResult = SimulationResult(integrator_tube,integrator_plate);

#   # Solve

#   ### Run the simulation and store data


@showprogress for t in tspan[1]:tstep:tspan[2]

    timemarching!(integrator_tube,integrator_plate,tstep)

    if (mod(integrator_plate.t,dt_record) < 1e-6) || (mod(-integrator_plate.t,dt_record) < 1e-6)
        store!(SimuResult,integrator_tube,integrator_plate)
    end

end

#   # Store data

save_path = "../numedata/solution.jld2"
save(save_path,"SimulationResult",SimuResult)

# ### take a peek at the solution (more at the PostProcessing notebook)

@gif for i in eachindex(SimuResult.tube_hist_t)
    plot(OHPTemp(),i,SimuResult,clim=(291.2,294.0))
end