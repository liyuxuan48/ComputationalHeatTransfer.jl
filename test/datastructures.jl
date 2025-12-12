using Statistics

"""
    Xtovec(Xp,dXdt) -> Vector

Transform vectors `Xp`, `dXdt` with the coordinates and velocities of the liquid slug interfaces into
a single state vector 
"""
function Xtovec(Xp::Array{Tuple{Float64,Float64},1},dXdt::Array{Tuple{Float64,Float64},1})
        
    Np = length(Xp)
    u = zeros(4*Np)
    
    for i = 1:Np
        # input Xp
        u[2*i-1] = Xp[i][1]
        u[2*i] = Xp[i][end]
        # input dXdt
        u[2*Np + 2*i-1] = dXdt[i][1]
        u[2*Np + 2*i] = dXdt[i][end]
    end
    
    return u
end


"""
    randomXp(L::Real,Lmin::Real,closedornot::Bool;numofslugs=DEFAULT_SLUGNUM,
                                                       chargeratio=DEFAULT_LIQUID_CHARGE_RATIO,
                                                       σ_charge=DEFAULT_SIGMA_CHARGE)->  X0,dXdt0,liquid_realratio

Transform vectors `Xp`, `dXdt` with the coordinates and velocities of the liquid slug interfaces into
a single state vector 
"""
DEFAULT_SLUGNUM = 30
DEFAULT_SIGMA_CHARGE = 0.01
DEFAULT_LIQUID_CHARGE_RATIO = 0.46
function randomXp(L::Float64,Lmin::Float64,closedornot::Bool;numofslugs=DEFAULT_SLUGNUM,
                                                       chargeratio=DEFAULT_LIQUID_CHARGE_RATIO,
                                                       σ_charge=DEFAULT_SIGMA_CHARGE)


    σ_persection = σ_charge*L/sqrt(numofslugs)

    L_perslug=L/numofslugs*chargeratio
    L_persection=L/numofslugs

    Ls = abs.((rand(numofslugs) .- 0.5).*sqrt(12).*σ_persection .+ L_perslug)

    Xp1s = zeros(numofslugs);
    Xp2s = deepcopy(Xp1s);

    if minimum(Ls) > Lmin && maximum(Ls) < L_persection
        if closedornot == true

            for i in eachindex(Xp1s)
                Xp1s[i] = (i-1)*L_persection
                Xp2s[i] = Xp1s[i] + Ls[i]
            end

            displacement = L*rand()

            Xp1s = mod.(Xp1s.+displacement,L)
            Xp2s = mod.(Xp2s.+displacement,L)

        # add openloop(starting from 0 for simplicity for now)
        elseif closedornot == false && numofslugs != 1
            for i in eachindex(Xp1s[1:end-1])
                Xp1s[i] = (i-1)*L_persection
                Xp2s[i] = Xp1s[i] + Ls[i]
            end
            Xp2s[end] = L
            Xp1s[end] = Xp2s[end] - Ls[end]

            displacement = 0.0
        end

    else
        error("Generation of random slugs failed")
    end

    X0 = map(tuple,Xp1s,Xp2s)
    dXdt0 = [zero.(X) for X in X0]
    liquid_realratio = sum(Ls)/L

    X0,dXdt0,liquid_realratio
end

@testset "Liquid and vapor arrays" begin

    L, Lmin, ϕ0 = 4.0 + rand(), 0.001, 0.6 + 0.1*rand()
    for closedornot in [true]
        X, dXdt, liquid_realratio = randomXp(L,Lmin,closedornot;chargeratio=ϕ0)

        # if this is a closed case, then shift it so that the last
        # slug crosses the line. Ensures more robust testing
        xshift = closedornot ? L - mean(X[end]) : 0.0
        X .= map(u -> (mod(u[1]+xshift,L),mod(u[2]+xshift,L)),X)        

        lslug = mod.(map(u -> u[2],X) - map(u -> u[1],X),L)
        @test all(lslug .> 0)

        lslug2 = ComputationalHeatTransfer.XptoLliquidslug(X,L)
        @test lslug == lslug2

        ϕ = sum(lslug)/L
        @test abs(ϕ - ϕ0) < 0.1
        @test ϕ ≈ liquid_realratio

        Xvapor = ComputationalHeatTransfer.getXpvapor(X,L,closedornot)
        lvap = mod.(map(u -> u[2],Xvapor) - map(u -> u[1],Xvapor),L)
        @test all(lvap .> 0)

        lvap2 = ComputationalHeatTransfer.XptoLvaporplug(X,L,closedornot)
        
        @test lvap == lvap2

        # Continuity test
        @test sum(lvap) + sum(lslug) ≈ L

        # Test checking positions in liquid slugs.
        np = length(X)

        # Liquid point should return true
        ir = rand(1:np-1)
        @test ComputationalHeatTransfer.ifamong(mean(X[ir]),X,L)

        # Vapor point should return false
        ir = rand(2:np)
        Xi, Xf = X[ir-1][2], X[ir][1]
        @test !ComputationalHeatTransfer.ifamong(0.5*(Xi+Xf),X,L)

        # Test creation of liquid slug arrays
        N = rand(200:300)
        θi = rand()
        Xarray, θarray = ComputationalHeatTransfer.constructXarrays(X,N,θi,L)
        # check that all points lie between 0 and L
        @test all(map(u -> all(0 .<= u .<= L),Xarray))
        # both arrays for a slug are the same length
        @test all(map((u,v) -> length(u)==length(v),Xarray,θarray))


        ## test assembly and disassembly 
        
        # test Xtovec, which was not in my original code
        dXdt .= [(rand(),rand()) for j in 1:np]
        u = Xtovec(X,dXdt)

        ir = rand(1:np)
        Xi, Xf = X[ir]
        dXi, dXf = dXdt[ir]
        @test u[2*ir] == Xf && u[2*ir-1] == Xi
        @test u[2*np + 2*ir] == dXf && u[2*np + 2*ir - 1] == dXi

        M = rand(np)
        δstart = rand(np)
        δend = rand(np)
        Lfilm_start = rand(np)
        Lfilm_end = rand(np)
        
        # Test that assembly and disassembly are inverses of each other
        u = ComputationalHeatTransfer.XMδLtovec(X,dXdt,M,δstart,δend,Lfilm_start,Lfilm_end)

        X2,dXdt2,M2,δstart2,δend2,Lfilm_start2,Lfilm_end2 = ComputationalHeatTransfer.vectoXMδL(u)

        @test X2 == X
        @test dXdt2 == dXdt
        @test M2 == M
        @test δstart2 == δstart
        @test δend2 == δend
        @test Lfilm_start2 == Lfilm_start
        @test Lfilm_end2 == Lfilm_end

        # # test a case in which there is one more X and dXdt element than the others
        # # This occurs in an open tube, e.g., vapor plugs at ends and one slug between them
        # # so that Xp contains (0,Xp1) and (Xpn,L) to describe the end vapor plugs.
        # X = [(rand(),rand()) for j in 1:np+1]
        # dXdt = [(rand(),rand()) for j in 1:np+1]
        # u = ComputationalHeatTransfer.XMδLtovec(X,dXdt,M,δstart,δend,Lfilm_start,Lfilm_end)
        # X2,dXdt2,M2,δstart2,δend2,Lfilm_start2,Lfilm_end2 = ComputationalHeatTransfer.vectoXMδL(u)

        # @test X2 == X
        # @test dXdt2 == dXdt
        # @test M2 == M
        # @test δstart2 == δstart
        # @test δend2 == δend
        # @test Lfilm_start2 == Lfilm_start
        # @test Lfilm_end2 == Lfilm_end

        

    end

end


#   # ASETS-II cases simulation

#   This is an example of a simulation package for conjugate heat transfer of an
#   oscillating heat pipe. SI units are used and units are emitted

#   ### What do we need to solve an OHP problem?

# 
#   **specify properties** : Solid property, Fluid property
# 
#   **set the geometries** : Computational domain, Heaters/Condensers, OHP shapes
# 
#   **construct the systems** : Fluid system(1D), HeatConduction system(2D)
# 
#   **initialize** : initialize the integrators and the data structs for saving
# 
#   **solve** : time marching to solve the two weakly coupled integrators
#   alternately
# 
#   **save/examine** : save the data for post-processing

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

#   ### Set up the evaporators and condensers

#   Right now, the OHPtype looks up a preset dictionary of OHP evaporators and
#   condensers.

#   You can also customize them in the OHP DIY notebook

OHPtype = "ASETS-II OHP 2 SMALL HEATER"
power = 10 + rand() # total heater power in watts
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


sys_plate = HeatConduction(params,Δx,xlim,ylim,Δt_max,qline=ohpgeom,qflux=eparams,qmodel=cparams)

@testset "Plate" begin
    ohp = sys_plate.qline
    
    @test ohp[1].body == ohpgeom.body                                    
end


@testset "Plate ADI solver" begin

    sys_plate_test = deepcopy(sys_plate)
    u_plate = newstate(sys_plate_test) .+ Tref # initialize plate T field to uniform Tref
    integrator_plate = init(u_plate,tspan,sys_plate_test) # construct integrator_plate

    ohp = sys_plate.qline
    
    @test ohp[1].body == ohpgeom.body

    #	Single eigenmode decay

    xylim = sys_plate_test.grid.xlim
    Δx = sys_plate_test.grid.Δx
    I0 = sys_plate_test.grid.I0

    xgrid = zero(temperature(integrator_plate))
    ygrid = zero(xgrid)
    for i in eachindex(xgrid[:,1])
        xgrid[i,:] .= xylim[1][1] .+ Δx*(i-1)
    end
    for i in eachindex(ygrid[1,:])
        ygrid[:,i] .= xylim[2][1] .+ Δx*(i-1)
    end

    temperature(integrator_plate) .= cos.(π .* xgrid ./ (xylim[1][2] .+ 0.5Δx)) .* cos.(2π .* ygrid ./ (xylim[2][2] .+ 0.5Δx))
    α =integrator_plate.p.params.α
    ρₛ = integrator_plate.p.params.ρ
    cₛ = integrator_plate.p.params.c

    T_hist = []
    Tgrid = temperature(integrator_plate);
    Tgrid_ini = deepcopy(Tgrid)
    tstep = 5e-3
    tend = 1.0
    ts = tstep:tstep:tend

    for i in ts
        push!(T_hist,deepcopy(ADI_newT!(Tgrid,sys_plate_test,tstep)))
    end


    j = rand(1:size(T_hist,1))
    Tratio_max = maximum(T_hist[j] ./ Tgrid_ini)
    Tratio_min = minimum(T_hist[j] ./ Tgrid_ini)

    # verify they all decay at the same rate
    @test isapprox(Tratio_max,Tratio_min,rtol=1e-10)

    Lx = xylim[1][2] .+ 0.5Δx
    Ly = xylim[2][2] .+ 0.5Δx
    
    # then verify maximum magnitude decay exponentially
    @test isapprox(log.(maximum.(T_hist) ./ maximum(Tgrid_ini)),-α * ((π/Lx)^2 + (2π/Ly)^2) .* ts,rtol=2e-3)

    # small heater case for one step
    T0 = Tref + 0.0 # uniform initial temperature

    temperature(integrator_plate) .= T0 # reset to uniform Tref
    T_hist = []
    Tgrid = temperature(integrator_plate);
    Tgrid_ini = deepcopy(Tgrid)
    tstep = 1e-2
    ts = tstep:tstep:tstep

    for i in ts
        push!(T_hist,deepcopy(ADI_timemarching!(Tgrid,sys_plate_test,tstep)))
    end

    # test mean temperature (energy conservation)

    T_hist_mean = mean(T_hist[1])
    ΔT_mean = T_hist_mean - T0
    ΔT_mean_analytical = power * tstep / (ρₛ * cₛ * 2Lx * 2Ly * plate_d)

    @test isapprox(ΔT_mean,ΔT_mean_analytical,atol=1e-12,rtol=1e-4)

    # test T rise/drop for one step with heater and condenser
    T0 = Tref + 10.0 # uniform initial temperature

    temperature(integrator_plate) .= T0 # reset to uniform Tref
    T_hist = []
    Tgrid = temperature(integrator_plate);
    Tgrid_ini = deepcopy(Tgrid)
    tstep = 1e-2
    ts = tstep:tstep:tstep

    for i in ts
        push!(T_hist,deepcopy(ADI_timemarching!(Tgrid,sys_plate_test,tstep)))
    end

    maximum_T = maximum(T_hist[1])
    minimum_T = minimum(T_hist[1])
    ΔT_max = maximum_T - T0
    ΔT_min = minimum_T - T0
    total_heater_area = 0.5inches*0.5inches;
    qe = power/total_heater_area
    ΔT_max_analytical = qe * tstep / (ρₛ * cₛ * plate_d)


    hc = sys_plate_test.qhdT[2].hc
    qc = -hc*(T0-Tref)
    ΔT_min_analytical = qc * tstep / (ρₛ * cₛ * plate_d)

    @test isapprox(ΔT_max,ΔT_max_analytical,atol=1e-12,rtol=5e-3)
    @test isapprox(ΔT_min,ΔT_min_analytical,atol=1e-12,rtol=5e-3)


    # test energy rise/drop for one step with small heater and line source

    sys_plate_test2 = ComputationalHeatTransfer.HeatConduction(params,Δx,xlim,ylim,Δt_max,qline=ohpgeom)
    # sys_tube = initialize_ohpsys(sys_plate_test2,p_fluid,power)
    u_plate2 = newstate(sys_plate_test2) .+ Tref # initialize plate T field to uniform Tref
    integrator_plate2 = init(u_plate2,tspan,sys_plate_test2) #

    xylim = sys_plate_test2.grid.xlim
    Lx = xylim[1][2] .+ 0.5Δx
    Ly = xylim[2][2] .+ 0.5Δx

    T0 = Tref + 0.0 # uniform initial temperature
    q1D_value = -2.0 # W/m line heat flux
    q1D = zeros(length(sys_plate_test2.qline[1].arccoord)) .+ q1D_value # W/m line heat flux
    # q1D = deepcopy(zero(sys_plate_test2.qline[1].arccoord) .+ 1.0) # W/m line heat flux
    L = sys_plate_test2.qline[1].arccoord[end]
    set_linesource_strength!(sys_plate_test2,q1D)

    temperature(integrator_plate2) .= T0 # reset to uniform Tref
    T_hist = []
    Tgrid = temperature(integrator_plate2);
    Tgrid_ini = deepcopy(Tgrid)
    tstep = 1e-2
    ts = tstep:tstep:tstep

    for i in ts
        push!(T_hist,deepcopy(ADI_timemarching!(Tgrid,sys_plate_test2,tstep)))
    end

    T_hist_mean = mean(T_hist[1])
    ΔT_mean = T_hist_mean - T0
    ΔT_mean_analytical = (-L*q1D_value) * tstep / (ρₛ * cₛ * 2Lx * 2Ly * plate_d)

    @test isapprox(ΔT_mean,ΔT_mean_analytical,atol=1e-12,rtol=1e-4)





    
    # temperature(sys_plate) .= 

                                                 
end

#   ### Create OHP inner channel system

#   sys_tube: fluid module system

sys_tube = initialize_ohpsys(sys_plate,p_fluid,power)

@testset "Film areas" begin

    d = sys_tube.tube.d
    δstart = sys_tube.vapor.δstart
    δend = sys_tube.vapor.δend
    dXdt = sys_tube.liquid.dXdt
    Ac = sys_tube.tube.Ac

    δdep = 0.05*rand()*d

    δarea_start = Ac .* (1 .- ((d .- 2*δstart) ./ d) .^ 2)
    δarea_end = Ac .* (1 .- ((d .- 2*δend) ./ d) .^ 2)
    δarea_dep = Ac .* (1 .- ((d .- 2*δdep) ./ d) .^ 2)
    @test all(δarea_start .== getδarea.(Ac,d,δstart))
    @test all(δarea_end .== getδarea.(Ac,d,δend))

    # No slugs are moving. Should only be equal to existing films
    Adep = getAdeposit(sys_tube,δdep)
    @test all(map((u,v) -> u[1]==v,Adep,δarea_end))
    @test all(map((u,v) -> u[2]==v,Adep,circshift(δarea_start,-1)))

    # Left slug interfaces are advancing. Should be equal to deposited film
    dXdt .= [(0.1,0.0) for i in eachindex(dXdt)]
    Adep = getAdeposit(sys_tube,δdep)
    @test all(map(u -> u[1]==δarea_dep,Adep))
    @test all(map((u,v) -> u[2]==v,Adep,circshift(δarea_start,-1)))

    # Right slug interfaces are advancing. Should be equal to deposited film
    dXdt .= [(0.0,-0.1) for i in eachindex(dXdt)]
    Adep = getAdeposit(sys_tube,δdep)
    @test all(map((u,v) -> u[1]==v,Adep,δarea_end))
    @test all(map(u -> u[2]==δarea_dep,Adep))

    vol = ComputationalHeatTransfer.getVolumevapor(sys_tube)
    @test all(vol .> 0)

    ρv = sys_tube.tube.PtoD.(sys_tube.vapor.P)
    @test ComputationalHeatTransfer.getMvapor(sys_tube) ≈ ρv.*vol

end

@testset "Mass functions" begin

    numofslugs = length(sys_tube.liquid.Xp)

    sys_tube.vapor.δstart .= zeros(numofslugs) .+ 1e-5 .+ 1e-5 .* rand(numofslugs) # initial velocity of the slugs
    sys_tube.vapor.δend .= zeros(numofslugs) .+ 1e-5 .+ 1e-5 .* rand(numofslugs) # initial velocity of the slugs
    

    d = sys_tube.tube.d
    δstart = sys_tube.vapor.δstart
    δend = sys_tube.vapor.δend
    Lfilm_start = sys_tube.vapor.Lfilm_start
    Lfilm_end = sys_tube.vapor.Lfilm_end
    Xp = sys_tube.liquid.Xp
    dXdt = sys_tube.liquid.dXdt
    Ac = sys_tube.tube.Ac
    L = sys_tube.tube.L
    ρₗ = sys_tube.liquid.ρ
    closedornot = sys_tube.tube.closedornot

    Lvaporplug = XptoLvaporplug(Xp,L,closedornot)
    Lliquidslug = XptoLliquidslug(Xp,L)
    Astart = getδarea(Ac,d,δstart)
    Aend = getδarea(Ac,d,δend)

    ρv = sys_tube.tube.PtoD.(sys_tube.vapor.P)

    # mass of vapor
    vol_vapor_analytical = Lvaporplug .* Ac .- Lfilm_start .* Astart .- Lfilm_end .* Aend
    M_vapor_analytical = ρv .* vol_vapor_analytical

    @test ComputationalHeatTransfer.getMvapor(sys_tube) ≈ M_vapor_analytical

    # mass of liquid
    vol_liquid_analytical = Lliquidslug .* Ac
    M_liquid_analytical = ρₗ .* vol_liquid_analytical

    @test ComputationalHeatTransfer.getMliquid(sys_tube) ≈ M_liquid_analytical

    # mass of films
    vol_film_start_analytical = Lfilm_start .* Astart
    M_film_start_analytical = ρₗ .* vol_film_start_analytical
    vol_film_end_analytical = Lfilm_end .* Aend
    M_film_end_analytical = ρₗ .* vol_film_end_analytical

    Mfilm_start,Mfilm_end = ComputationalHeatTransfer.getMfilm(sys_tube)

    @test Mfilm_start ≈ M_film_start_analytical
    @test Mfilm_end ≈ M_film_end_analytical

end

# some threshold values need to be changed for other applications
@testset "Hfilm" begin

    δmin = sys_tube.vapor.δmin;
    δthreshold = 5e-6
    δmax = 1e-4

    kₗ = sys_tube.vapor.k
    Hᵥ = sys_tube.vapor.Hᵥ
    δs = [1e-6;3e-6;1e-5;1.5e-4]

    @test ComputationalHeatTransfer.Hfilm.(δs,[sys_tube]) ≈ [0.0;(δs[2]-δmin)*(kₗ/δthreshold - Hᵥ)/(δthreshold-δmin);kₗ/δs[3];0.5*kₗ/δmax+1e-6]
end