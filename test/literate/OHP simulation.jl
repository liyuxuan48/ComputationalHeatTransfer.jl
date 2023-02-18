# # ASETS-II

# This is an example of a simulation package for conjugate heat transfer of an oscillating heat pipe. SI units are used and units are emitted

# ## Packages
# If you need to update the ComputationalHeatTransfer package, I suggest you activate the package's path as your current environment. If you don't need to modify the code within the packge, you don't have to run this block.


using Pkg
Pkg.activate(dirname(dirname(pwd()))) # using current environment for development

# Firstly, let's import the necessary packages, you may need to install them for the first time.


using ComputationalHeatTransfer # our main package
using LaTeXStrings # for latex strings
using JLD2 # for file input/output
using ProgressMeter # to have a progessbar when runing the simulation
using XLSX # for reading experimental data in CSV format
using Plots # for plotting
gr()  #ploting backend (the fastest one)


# ## Control Console

# This block contains the parameters I reckon could be tuned to match the experimental data. 
# As I am still tuning, they are placed here for convenience.


OHPtype = "ASETS-II OHP 2 LARGE HEATER"
power = 40 # total heater power in watts

hc = 3000.0 #condenser heat transfer coefficient
Rn = 3e-6 # nucleation site radius
δfilm = 2e-5 # initial film thickness
ad_fac = 1.3 # film thickness factor
plate_d = 1.5e-3; # plate thickness
η₊ = 0.2 # η+ for the film model
η₋ = 0.15 # η- for the film model(sugguested smaller than η+)

# # Properies

# ## Solid Physical parameters

# params is the HeatConductionParameters for the plate material (aluminum).


ρₛ = 2730; # density
cₛ  = 8.93e02; # specific heat
kₛ  = 1.93e02; # heat conductivity
params = HeatConductionParameters(ρₛ ,cₛ ,kₛ ,thickness=plate_d)


# ## Fluid Physical parameters

# p_fluid contains the vapor and liquid properties at a constant reference temperature. Noted that the vapor pressure and the vapor density will be functions of temperatures during the simulation, other properties are extracted from p_fluid as an approximate value.

Tref = 291.2 # reference temperature
fluid_type = "Butane"
p_fluid = SaturationFluidProperty(fluid_type,Tref)

# # Plate Conduction Part

# ### Geometry parameters

# The 2D domain is of rectangular shape (slightly different from ASETS-II). In the future it can be of arbitrary shape using the immersedlayers.jl package.


Lx = 0.1524; # plate size x [m]
Ly = 0.0648; # plate size y [m]
xlim = (-Lx/2,Lx/2) # plate x limits
ylim = (-Ly/2,Ly/2) # plate y limits

# ### Set mesh size and maximum time step for plate heat conduction

# Δx is controlled by Δx = α*gridPe and set having the same order of magitute of tube diameter 1e-3. Fourier number is used to give a safety "cap" of time step you can choose in the fluid module


Δx,Δt_max = setstepsizes(params.α,gridPe=8.0,fourier=0.3) 

# ### Set up the evaporators and condensers
# In one sided condenser case, the "adiabatic" side is not completely adiabatic in the code. Instead it is a fraction "hc2ratio" of the regular condenser htc. And it is a tunable parameter for now as it is a representaion of the insulation material:)


Tc = Tref; # condenser temperature
eparams,cparams = OHPConfiguration(OHPtype,power,Tc,hc,Δx,hc2ratio=1/50);

# ### Set up OHP channel's shape
# construct_ohp_curve is a built-in function that generates two arrays: x that contains all x values of the discrete points, and y contains all y values. x and y have the same length. You can also customize this function to generate an OHP shape of your choice as long as they produce x and y.


x, y = construct_ohp_curve("ASETS",Δx) # get x and y coordinates for the channel
ohp = BasicBody(x,y) # build a BasicBody based on x,y
ohpgeom = ComputationalHeatTransfer.LineSourceParams(ohp) # build a line heat source based on BasicBody

# This is a exmaple of the compuational domain (the box) and the OHP channel serpentine (in blue)

plot(ohp,fillalpha=0,linecolor=:blue,xlims=xlim,ylims=ylim,framestyle = :box)

# ### Create HeatConduction system
# The solid module dealing with the 2D conduction, evaporator, condenser, and the OHP line heat source is constructed here.


sys = HeatConduction(params,Δx,xlim,ylim,Δt_max,qline=ohpgeom,qflux=eparams,qmodel=cparams)

# ### Create OHP inner channel system

# sys_tube: fluid module system

# u       : fluid module state vector

# cbst    : fluid module callback function


sys_tube,u,cbst = initialize_ohpsys(OHPtype,fluid_type,sys,p_fluid,Tref,δfilm,η₊,η₋,Rn,ad_fac);

# ### set time step


tspan = (0.0, 1.0); # start time and end time
dt_record = 0.2     # saving time interval
num_data = (tspan[2] - tspan[1]) / dt_record # calculate the number of saving data points

global tstep = 5e-4        # actrual time marching step

# ### combine inner tube and plate together


u_plate = newstate(sys) .+ Tref # initialize T field
integrator_plate = init(u_plate,tspan,sys) # construct integrator_plate


prob = ODEProblem(ODE_innertube, u, tspan, sys_tube) # construct integrator_tube problem
integrator_tube = init(prob, RK4(),save_on=false, dt=tstep, callback=cbst,maxiters=1e10); # construct integrator_tube

# ### initialize arrays for saving


boil_hist=[]
plate_T_hist = []
tube_hist_u  = []
tube_hist_t = []
tube_hist_θwall = []

# ### Run the simulation and store data


@showprogress for t in tspan[1]:tstep:tspan[2]

    timemarching!(integrator_tube,integrator_plate,tstep)
    
    if (mod(integrator_plate.t,dt_record) < 1e-6) || (mod(-integrator_plate.t,dt_record) < 1e-6)
        push!(plate_T_hist,deepcopy(temperature(integrator_plate))); 
        push!(tube_hist_θwall,deepcopy(integrator_tube.p.wall.θarray))
        push!(tube_hist_u,deepcopy(integrator_tube.u));
        push!(tube_hist_t,deepcopy(integrator_tube.t));
        integrator_plate_temp = deepcopy(integrator_plate)
        println(integrator_plate.t)
    end
    
end
# integrator_plate.t

# ### Convert state vector to system, for post-processing conveience


sysfinal = []
for i = 1:length(tube_hist_u)
    sysfinaltemp = deepcopy(getcurrentsys(tube_hist_u[i],integrator_tube.p))
    sysfinaltemp.wall.θarray = tube_hist_θwall[i]
    push!(sysfinal, sysfinaltemp)
end

# ## Store data


boil_station_hist = [elem[1] for elem in boil_hist]
boil_t_hist = [elem[2] for elem in boil_hist]
boil_data = Array(reshape([boil_station_hist;boil_t_hist], length(boil_t_hist),2));


save(dirname(dirname(dirname(pwd())))*"/OHPnume/OHP2_40W_large.jld2",
    "tube_hist_θwall", tube_hist_θwall,"tube_hist_u", tube_hist_u,"tube_hist_t", tube_hist_t,
    "integrator_tube",integrator_tube,"plate_T_hist", plate_T_hist,
    "integrator_plate", integrator_plate, "boil_data", boil_data,"ohp",ohp)
