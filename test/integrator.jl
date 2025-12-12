using Statistics

heaviside(x::AbstractFloat) = ifelse(x < 0, zero(x), ifelse(x > 0, one(x), oftype(x,0.0)))


using ComputationalHeatTransfer # our main package

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

#   ### Set up the evaporators and condensers

#   Right now, the OHPtype looks up a preset dictionary of OHP evaporators and
#   condensers.

#   You can also customize them in the OHP DIY notebook

OHPtype = "ASETS-II OHP 2 LARGE HEATER"
power = 40 # total heater power in watts
Tc = Tref; # condenser temperature
eparams,cparams = OHPConfiguration(OHPtype,power,Tc,Δx);

#   ### Set up OHP channel's shape

#   constructohpcurve is a built-in function that generates two arrays: x that
#   contains all x values of the discrete points, and y contains all y values. x
#   and y have the same length.

#   You can also customize this function to generate an OHP shape of your choice
#   as long as they produce x array and y array.

x, y = construct_ohp_curve("ASETS",Δx) # get x and y coordinates for the channel
ohp = BasicBody(x,y) # build a BasicBody based on x,y

ohpgeom = ComputationalHeatTransfer.LineSourceParams(ohp) # build a line heat source based on BasicBody

#   # Construct the systems

#   ### Create HeatConduction system

#   The solid module dealing with the 2D conduction, evaporator, condenser, and
#   the OHP line heat source is constructed here.

sys_plate = HeatConduction(params,Δx,xlim,ylim,Δt_max,qline=ohpgeom,qflux=eparams,qmodel=cparams)

#   ### Create OHP inner channel system

#   sys_tube: fluid module system

# test 1, boiling callbacks
@testset "merging callbacks" begin

    tspan = (0.0, 5.0); # start time and end time
    dt_record = 0.01   # saving time interval

    tstep = 1e-3     # actrual time marching step


    A_SMALL_FRAC = 1e-2

    sys_tube = initialize_ohpsys(sys_plate,p_fluid,power)
    u_tube_1 = newstate(sys_tube) # initialize OHP tube 
    integrator_tube_1 = init(u_tube_1,tspan,deepcopy(sys_tube)); # construct integrator_tube
    @test merging_condition(integrator_tube_1.u,integrator_tube_1.t,integrator_tube_1) == false

    L_newbubble= sys_tube.wall.L_newbubble
    L = sys_tube.tube.L

    Xp_old = sys_tube.liquid.Xp
    Xp_new = deepcopy(Xp_old)

    dXp = mod(Xp_new[2][1] - Xp_new[1][end],L) - 0.4*L_newbubble
    @test  dXp > 0.0

    Xp_new[2] =  (Xp_new[2][1] - dXp,Xp_new[2][2] - dXp) # make two vapor plugs close enough to merge

    Lfilm_start= sys_tube.vapor.Lfilm_start
    Lfilm_end= sys_tube.vapor.Lfilm_end

    Lfilm_start[2] = A_SMALL_FRAC*L_newbubble
    Lfilm_end[2] = A_SMALL_FRAC*L_newbubble

    sys_tube.liquid.Xp = deepcopy(Xp_new)
    sys_tube.vapor.Lfilm_start = deepcopy(Lfilm_start)
    sys_tube.vapor.Lfilm_end = deepcopy(Lfilm_end)

    u_tube_2 = newstate(sys_tube) # initialize OHP tube 
    integrator_tube_2 = init(u_tube_2,tspan,deepcopy(sys_tube)); # construct integrator_tube
    @test merging_condition(integrator_tube_2.u,integrator_tube_2.t,integrator_tube_2) == true

    # println("Before modifying Xp: ", sys_tube.liquid.Xp)
    # println("After modifying Xp: ", Xp_new)

    #   ### combine inner tube and plate together

    u_plate = newstate(sys_plate) .+ Tref # initialize plate T field to uniform Tref
    integrator_plate = init(u_plate,tspan,sys_plate) # construct integrator_plate

    u_tube = newstate(sys_tube) # initialize OHP tube 
    integrator_tube = init(u_tube,tspan,sys_tube); # construct integrator_tube

    SimuResult = SimulationResult(integrator_tube,integrator_plate);


    # println(integrator_tube.p.liquid.Xp)


    p_old = deepcopy(integrator_tube.p)

    Mvapor_old = sum(getMvapor(p_old))
    Mfilm_old = sum(sum.(getMfilm(p_old)))
    Mliquid_old = sum(getMliquid(p_old))

    merging_affect!(integrator_tube)
    getcurrentsys!(integrator_tube.u,integrator_tube.p)

    p_new = deepcopy(integrator_tube.p)

    @test length(p_new.liquid.Xp) == length(p_old.liquid.Xp) - 1
    @test isapprox(sum(XptoLliquidslug(p_new.liquid.Xp,p_new.tube.L)), sum(XptoLliquidslug(p_old.liquid.Xp,p_old.tube.L)), rtol=2e-3)
    @test p_new.liquid.Xp[2:end] == p_old.liquid.Xp[3:end]
    @test p_new.liquid.dXdt[2:end] == p_old.liquid.dXdt[3:end]


    L_liquidslug_old = XptoLliquidslug(p_old.liquid.Xp,p_old.tube.L)
    L_liquidslug_new = XptoLliquidslug(p_new.liquid.Xp,p_new.tube.L)
    @test p_new.liquid.dXdt[1][1] ≈ (p_old.liquid.dXdt[1][1] * L_liquidslug_old[1] + p_old.liquid.dXdt[2][1] * L_liquidslug_old[2])/(L_liquidslug_old[1]+L_liquidslug_old[2]) 

    Mvapor_new = sum(getMvapor(p_new))
    Mfilm_new = sum(sum.(getMfilm(p_new)))
    Mliquid_new = sum(getMliquid(p_new))

    @test isapprox(Mvapor_new+Mfilm_new+Mliquid_new, Mvapor_old+Mfilm_old+Mliquid_old, rtol=1e-10)



    # @showprogress for t in tspan[1]:tstep:tspan[1]

    #     timemarching!(integrator_tube,integrator_plate,tstep)

    #     if (mod(integrator_plate.t,dt_record) < 1e-6) || (mod(-integrator_plate.t,dt_record) < 1e-6)
    #         store!(SimuResult,integrator_tube,integrator_plate)
    #     end

    # end
end

@testset "boiling callbacks" begin

    numofslugs = 2

    sys_tube = initialize_ohpsys(sys_plate,p_fluid,power,slugnum=numofslugs,nucleatenum=2,
        boil_waiting_time=1e-2)

    Xp = sys_tube.liquid.Xp
    L = sys_tube.tube.L

    Xp1_mid = mod(Xp[1][1] + mod((Xp[1][2]- Xp[1][1]),L)/2,L)
    Xp2_vapor_mid = mod(Xp[1][2] + mod((Xp[2][1]- Xp[1][2]),L)/2,L)

    sys_tube.wall.Xstations=[Xp1_mid,Xp2_vapor_mid]


    tspan = (0.0, 1.0); # start time and end time
    dt_record = 1.0   # saving time interval

    tstep = 1e-2     # actrual time marching step

    Rn = sys_tube.wall.Rn
    d = sys_tube.tube.d
    TtoP = sys_tube.tube.TtoP
    ΔT = RntoΔT(Rn,Tref,fluid_type,d,TtoP)
    

    u_plate = newstate(sys_plate) .+ Tref .+ 10.0# initialize plate T field to uniform Tref
    integrator_plate = init(u_plate,tspan,sys_plate) # construct integrator_plate

    u_tube = newstate(sys_tube) # initialize OHP tube 
    integrator_tube = init(u_tube,tspan,deepcopy(sys_tube)); # construct integrator_tube

    SimuResult = SimulationResult(integrator_tube,integrator_plate);

    timemarching!(integrator_tube,integrator_plate,tstep)

    @test integrator_tube.p.wall.boiltime_stations == [0.0,tstep]

    timemarching!(integrator_tube,integrator_plate,tstep)

    @test integrator_tube.p.wall.boiltime_stations == [2*tstep,2*tstep]


    u_plate = newstate(sys_plate) .+ Tref .+ 1.1ΔT# initialize plate T field to uniform Tref
    integrator_plate = init(u_plate,tspan,sys_plate) # construct integrator_plate


    u_tube = newstate(sys_tube) # initialize OHP tube 
    integrator_tube = init(u_tube,tspan,deepcopy(sys_tube)); # construct integrator_tube
    integrator_tube.p.wall.θarray = temperature_linesource(integrator_plate)
    

    SimuResult = SimulationResult(integrator_tube,integrator_plate);

    boiling_affect!(integrator_tube)
    getcurrentsys!(integrator_tube.u,integrator_tube.p)

    @test integrator_tube.p.liquid.Xp == sys_tube.liquid.Xp

    integrator_tube.t = integrator_tube.p.wall.boil_interval + tstep

    p_old = deepcopy(integrator_tube.p)
    Mvapor_old = sum(getMvapor(p_old))
    Mfilm_old = sum(sum.(getMfilm(p_old)))
    Mliquid_old = sum(getMliquid(p_old))

    boiling_affect!(integrator_tube)
    getcurrentsys!(integrator_tube.u,integrator_tube.p)

    @test length(integrator_tube.p.liquid.Xp) == length(sys_tube.liquid.Xp) + 1


    p_new = deepcopy(integrator_tube.p)
    Mvapor_new = sum(getMvapor(p_new))
    Mfilm_new = sum(sum.(getMfilm(p_new)))
    Mliquid_new = sum(getMliquid(p_new))
    @test isapprox(Mvapor_new+Mfilm_new+Mliquid_new, Mvapor_old+Mfilm_old+Mliquid_old, rtol=1e-4)
    @test p_old.liquid.dXdt[1][1] == p_new.liquid.dXdt[1][1] == p_new.liquid.dXdt[2][1]

    u_plate = newstate(sys_plate) .+ Tref .+ 0.9ΔT# initialize plate T field to uniform Tref
    integrator_plate = init(u_plate,tspan,sys_plate) # construct integrator_plate
    u_tube = newstate(sys_tube) # initialize OHP tube 
    integrator_tube = init(u_tube,tspan,deepcopy(sys_tube)); # construct integrator_tube
    integrator_tube.p.wall.θarray = temperature_linesource(integrator_plate)
    integrator_tube.t = integrator_tube.p.wall.boil_interval + tstep

    boiling_affect!(integrator_tube)
    getcurrentsys!(integrator_tube.u,integrator_tube.p)

    @test integrator_tube.p.liquid.Xp == sys_tube.liquid.Xp

end



@testset "fixdx callbacks" begin



    tspan = (0.0, 1.0); # start time and end time
    dt_record = 1.0   # saving time interval

    tstep = 1e-2     # actrual time marching step


 sys_tube = initialize_ohpsys(sys_plate,p_fluid,power)

u_tube = newstate(sys_tube) # initialize OHP tube 
integrator_tube = init(u_tube,tspan,deepcopy(sys_tube)); # construct integrator_tube

@test fixdx_condition(integrator_tube.u,integrator_tube.t,integrator_tube) == false

sys_tube.liquid.Xarrays[1] = ComputationalHeatTransfer.constructoneXarray(integrator_tube.p.liquid.Xp[1],3*length(integrator_tube.p.liquid.Xarrays[1]),integrator_tube.p.tube.L)
sys_tube.liquid.θarrays[1] = sys_tube.liquid.Xarrays[1] .* 0 .+ Tref
u_tube = newstate(sys_tube) # initialize OHP tube 
integrator_tube = init(u_tube,tspan,deepcopy(sys_tube)); # construct integrator_tube


@test fixdx_condition(integrator_tube.u,integrator_tube.t,integrator_tube) == true

fixdx_affect!(integrator_tube)
getcurrentsys!(integrator_tube.u,integrator_tube.p)

@test length(integrator_tube.p.liquid.Xarrays[1]) < length(sys_tube.liquid.Xarrays[1])
@test fixdx_condition(integrator_tube.u,integrator_tube.t,integrator_tube) == false

end

@testset "slugbc callbacks" begin



    tspan = (0.0, 1.0); # start time and end time
    dt_record = 1.0   # saving time interval

    tstep = 1e-2     # actrual time marching step


 sys_tube = initialize_ohpsys(sys_plate,p_fluid,power)
 
u_tube = newstate(sys_tube) # initialize OHP tube 
integrator_tube = init(u_tube,tspan,deepcopy(sys_tube)); # construct integrator_tube
# for i in eachindex(integrator_tube.p.vapor.P)
#     integrator_tube.p.liquid.θarrays[i] = integrator_tube.p.liquid.θarrays[i] .+ rand() .* 1e3
# end

integrator_tube.p.vapor.P = integrator_tube.p.vapor.P + rand(length(integrator_tube.p.vapor.P)) .* 1e3


@test slugbc_condition(integrator_tube.u,integrator_tube.t,integrator_tube) == true # always true


T_first_pt = [x[1] for x in integrator_tube.p.liquid.θarrays]
T_last_pt = [x[end] for x in integrator_tube.p.liquid.θarrays]

PtoT = integrator_tube.p.tube.PtoT

@test !isapprox(PtoT.(integrator_tube.p.vapor.P),T_first_pt, rtol=1e-10)
@test !isapprox(PtoT.(integrator_tube.p.vapor.P),circshift(T_last_pt,1), rtol=1e-10)

slugbc_affect!(integrator_tube)
getcurrentsys!(integrator_tube.u,integrator_tube.p)

@test isapprox(PtoT.(integrator_tube.p.vapor.P),T_first_pt, rtol=1e-10)
@test isapprox(PtoT.(integrator_tube.p.vapor.P),circshift(T_last_pt,1), rtol=1e-10)

end

# test 1, boiling callbacks
@testset "vapormerging callbacks" begin

    tspan = (0.0, 5.0); # start time and end time
    dt_record = 0.01   # saving time interval

    tstep = 1e-3     # actrual time marching step


    A_SMALL_FRAC = 1e-2

    sys_tube = initialize_ohpsys(sys_plate,p_fluid,power)
    u_tube_1 = newstate(sys_tube) # initialize OHP tube 
    integrator_tube_1 = init(u_tube_1,tspan,deepcopy(sys_tube)); # construct integrator_tube
    @test vaporMergingCondition(integrator_tube_1.u,integrator_tube_1.t,integrator_tube_1) == false

    L_newbubble= sys_tube.wall.L_newbubble
    L = sys_tube.tube.L

    Xp_old = sys_tube.liquid.Xp
    Xp_new = deepcopy(Xp_old)

    dXp = mod(Xp_new[2][2] - Xp_new[2][1],L) - 0.4*L_newbubble
    Xp_new[2] =  (Xp_new[2][1], Xp_new[2][2] - dXp) # make two liquid slugs close enough to merge

    sys_tube.liquid.Xp = deepcopy(Xp_new)

    u_tube_2 = newstate(sys_tube) # initialize OHP tube 
    integrator_tube_2 = init(u_tube_2,tspan,deepcopy(sys_tube)); # construct integrator_tube

    @test vaporMergingCondition(integrator_tube_2.u,integrator_tube_2.t,integrator_tube_2) == true

#     # println("Before modifying Xp: ", sys_tube.liquid.Xp)
#     # println("After modifying Xp: ", Xp_new)

#     #   ### combine inner tube and plate together

#     u_plate = newstate(sys_plate) .+ Tref # initialize plate T field to uniform Tref
#     integrator_plate = init(u_plate,tspan,sys_plate) # construct integrator_plate

    u_tube = newstate(sys_tube) # initialize OHP tube 
    integrator_tube = init(u_tube,tspan,sys_tube); # construct integrator_tube

#     SimuResult = SimulationResult(integrator_tube,integrator_plate);


#     # println(integrator_tube.p.liquid.Xp)


    p_old = deepcopy(integrator_tube.p)

    Mvapor_old = sum(getMvapor(p_old))
    Mfilm_old = sum(sum.(getMfilm(p_old)))
    Mliquid_old = sum(getMliquid(p_old))

    vaporMergingAffect!(integrator_tube)
    getcurrentsys!(integrator_tube.u,integrator_tube.p)

    p_new = deepcopy(integrator_tube.p)

    @test length(p_new.liquid.Xp) == length(p_old.liquid.Xp) - 1
    @test isapprox(sum(XptoLliquidslug(p_new.liquid.Xp,p_new.tube.L)), sum(XptoLliquidslug(p_old.liquid.Xp,p_old.tube.L)), rtol=2e-3)

    Mvapor_new = sum(getMvapor(p_new))
    Mfilm_new = sum(sum.(getMfilm(p_new)))
    Mliquid_new = sum(getMliquid(p_new))

    @test isapprox(Mvapor_new+Mfilm_new+Mliquid_new, Mvapor_old+Mfilm_old+Mliquid_old, rtol=1e-10)

end

@testset "weakly coupled time marching" begin

    sys_tube = initialize_ohpsys(sys_plate,p_fluid,power)

    tspan = (0.0, 1.0); # start time and end time
    dt_record = 1.0   # saving time interval

    tstep = 1e-3     # actrual time marching step

    u_plate = newstate(sys_plate) .+ Tref# initialize plate T field to uniform Tref
    integrator_plate = init(u_plate,tspan,sys_plate) # construct integrator_plate

    u_tube = newstate(sys_tube) # initialize OHP tube 
    integrator_tube = init(u_tube,tspan,deepcopy(sys_tube)); # construct integrator_tube

    SimuResult = SimulationResult(integrator_tube,integrator_plate);

    @test all(isapprox.(integrator_tube.p.wall.θarray,Tref, rtol=1e-10))

    timemarching!(integrator_tube,integrator_plate,tstep)

    @test all(isapprox.(integrator_tube.p.wall.θarray,Tref, rtol=1e-10))

    timemarching!(integrator_tube,integrator_plate,tstep)

    @test !all(isapprox.(integrator_tube.p.wall.θarray,Tref, rtol=1e-10))

    store!(SimuResult,integrator_tube,integrator_plate)

    # make sure after time marching, the system state is up-to-date with the latest u
    p_old = deepcopy(integrator_tube.p)
    systemp = ComputationalHeatTransfer.getcurrentsys!(integrator_tube.u,integrator_tube.p)

    @test all(p_old.vapor.P .== systemp.vapor.P)
    @test SimuResult.tube_hist_θwall[end] == integrator_tube.p.wall.θarray
    @test SimuResult.tube_hist_u[end] == integrator_tube.u
    @test SimuResult.tube_hist_t[end] == integrator_tube.t
end