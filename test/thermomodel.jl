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


# SimuResult = SimulationResult(integrator_tube,integrator_plate);

#   # Solve

#   ### Run the simulation and store data


# test 1, 30 slugs, no heat transfer, given an initial velocity
@testset "dynamicsmodel test1" begin
    numofslugs = 30
    sys_tube = initialize_ohpsys(sys_plate,p_fluid,power,slugnum=numofslugs,ηplus=rand())
    sys_tube.liquid.dXdt = [zero.(X) .+ 1.0 for X in sys_tube.liquid.dXdt]# initial velocity of the slugs
    sys_tube.vapor.δstart .= zeros(numofslugs) .+ 1e-5 # initial velocity of the slugs

    u_tube = newstate(sys_tube) # initialize OHP tube
    # integrator_tube = init(u_tube,tspan,sys_tube); # construct integrator_tube

    σ = sys_tube.liquid.σ
    L = sys_tube.tube.L
    ρₗ = sys_tube.liquid.ρ
    μₗ = sys_tube.liquid.μₗ
    ad_fac = sys_tube.vapor.ad_fac
    d = sys_tube.tube.d
    Ac = sys_tube.tube.Ac
    peri = sys_tube.tube.peri
    Xp = sys_tube.liquid.Xp
    dXdt = sys_tube.liquid.dXdt
    δend = sys_tube.vapor.δend
    # characteristic bulk velocities for each liquid slug
    V = [mean(elem) for elem in sys_tube.liquid.dXdt]
    # get a characteristic Capilarry number based on the average velocities
    Vavg = mean(abs.(V))
    Ca = getCa.(μₗ,σ,Vavg)

    δdep = Catoδ(d,Ca,adjust_factor=ad_fac)

    Lliquidslug = XptoLliquidslug(Xp,L)

    Adeposit = getAdeposit(sys_tube,δdep)
    Adeposit_left = [elem[1] for elem in Adeposit]
    Adeposit_right = [elem[2] for elem in Adeposit]
    
# The first test is to check the case where there is no heat transfer and an initial velocity
    uu_test1 = dynamicsmodel(u_tube[1:9*numofslugs],sys_tube)

    @test all(uu_test1[1:2:2*numofslugs-1] .≈ Ac ./ (Ac .- Adeposit_left) .*  V) # velocity should remain the same
    @test all(uu_test1[2:2:2*numofslugs]   .≈ Ac ./ (Ac .- Adeposit_right) .*  V) # velocity should remain the same

    # get differential equation factors
    lhs = ρₗ*Ac .* Lliquidslug

    # analytical solution for dXdt term (friction)
    Re = ρₗ .* abs.(V) .* d ./ μₗ
    f_coefficient = f_churchill.(Re)
    dXdt_to_stress = -1/8 .* f_coefficient .* ρₗ .* V .* abs.(V)
    rhs_dXdt = peri .* Lliquidslug .* dXdt_to_stress ./ lhs

    # analytical solution for dLdt term (mass conservation)
    dVdt_nodLdt = dXdt_to_stress*peri .* Lliquidslug ./ lhs # if dLdt = 0
    dLdt = (Ac ./ (Ac .- Adeposit_right) - Ac ./ (Ac .- Adeposit_left)) .*  V
    rhs_dLdt = -ρₗ .* Ac .* V .*  dLdt ./ lhs
    
    @test all(uu_test1[2*numofslugs+1:2:4*numofslugs-1] .≈ rhs_dXdt .+ rhs_dLdt)

    @test all(uu_test1[2*numofslugs+2:2:4*numofslugs]   .== uu_test1[2*numofslugs+1:2:4*numofslugs-1])

    @test all(isapprox.(uu_test1[4*numofslugs+1:5*numofslugs],0.0,atol=1e-12)) #dMdt = 0

    @test all(isapprox.(uu_test1[5*numofslugs+1:6*numofslugs],0.0,atol=1e-12)) #dδstart/dt = 0

    F_end = ρₗ .* getδarea.(Ac,d,δend)
    F2_end = ρₗ .* getδarea.(Ac,d,δdep)

    peri_end = peri .* (d .- 2δend)/d
    C_end = ρₗ .* peri_end
    Lfilm_end = sys_tube.vapor.Lfilm_end

    @test all(uu_test1[6*numofslugs+1:7*numofslugs] .≈ (F2_end .- F_end) .* Ac ./ (Ac .- Adeposit_left) .*  V ./ C_end ./ Lfilm_end)

    @test all(uu_test1[7*numofslugs+1:8*numofslugs] .≈ -uu_test1[2:2:2*numofslugs])

    @test all(uu_test1[8*numofslugs+1:9*numofslugs] .≈ uu_test1[1:2:2*numofslugs-1])
                                                 
end

# test 2, 2 slugs, with heat transfer, given initial temperature differences, no initial velocity
@testset "dynamicsmodel test2" begin

    ηplusvalue = rand()
    ηminusvalue = 0.0
    numofslugs = 2
    sys_tube = initialize_ohpsys(sys_plate,p_fluid,power,slugnum=numofslugs,ηplus=ηplusvalue)

    ΔT = rand()
    ΔT_array = [ΔT,-ΔT]
    sys_tube.vapor.P[1] = sys_tube.tube.TtoP(Tref-ΔT_array[1])
    sys_tube.vapor.P[2] = sys_tube.tube.TtoP(Tref-ΔT_array[2])

    closedornot = sys_tube.tube.closedornot
    P1 = sys_tube.vapor.P[1]
    P2 = sys_tube.vapor.P[2]
    L = sys_tube.tube.L
    Ac = sys_tube.tube.Ac
    d = sys_tube.tube.d
    ρₗ = sys_tube.liquid.ρ
    Xp = sys_tube.liquid.Xp
    Hfg = sys_tube.tube.PtoHfg.(sys_tube.vapor.P)
    k = sys_tube.vapor.k
    peri = sys_tube.tube.peri
    δstart = sys_tube.vapor.δstart
    δend = sys_tube.vapor.δend
    Lfilm_start = sys_tube.vapor.Lfilm_start
    Lfilm_end = sys_tube.vapor.Lfilm_end

    u_tube = newstate(sys_tube) # initialize OHP tube

    uu_test2 = dynamicsmodel(u_tube[1:9*numofslugs],sys_tube)


    @test all(isapprox.(uu_test2[1:2:2*numofslugs-1],0.0,atol=1e-12)) # velocity should remain the same
    @test all(isapprox.(uu_test2[2:2:2*numofslugs],0.0,atol=1e-12)) # velocity should remain the same

    Lliquidslug = XptoLliquidslug(Xp,L)
    lhs = ρₗ*Ac .* Lliquidslug
    @test all(uu_test2[2*numofslugs+1:2:4*numofslugs-1] .≈ [1,-1] .* (P1-P2)*Ac ./ lhs) # velocity should remain the same
    @test all(uu_test2[2*numofslugs+2:2:4*numofslugs]  .== uu_test2[2*numofslugs+1:2:4*numofslugs-1])

    # analytical solution for dMdt_latent
    dMdt_latent_start_analytical = Lfilm_start .* peri .* k ./ δstart ./ Hfg .* ΔT_array
    dMdt_latent_end_analytical = Lfilm_end .* peri .* k ./ δend ./ Hfg .* ΔT_array
    dMdt_latent_start_positive_analytical = heaviside.(dMdt_latent_start_analytical) .* dMdt_latent_start_analytical
    dMdt_latent_end_positive_analytical = heaviside.(dMdt_latent_end_analytical) .* dMdt_latent_end_analytical
    dMdt_latent_start_negative_analytical = heaviside.(-dMdt_latent_start_analytical) .* dMdt_latent_start_analytical
    dMdt_latent_end_negative_analytical = heaviside.(-dMdt_latent_end_analytical) .* dMdt_latent_end_analytical

    # run the dynamics model to get numerical dMdt_latent
    Xpvapor = getXpvapor(Xp,L,closedornot)
    dMdt_latent_start,dMdt_latent_end,dMdt_latent_start_positive,dMdt_latent_end_positive = dMdtdynamicsmodel(Xpvapor,sys_tube)
    dMdt_latent_start_negative = dMdt_latent_start .- dMdt_latent_start_positive
    dMdt_latent_end_negative = dMdt_latent_end .- dMdt_latent_end_positive

    @test all(isapprox.(dMdt_latent_start,dMdt_latent_start_analytical,rtol=4e-3)) 
    @test all(isapprox.(dMdt_latent_end,dMdt_latent_end_analytical,rtol=4e-3)) 
    @test all(isapprox.(dMdt_latent_start_positive,dMdt_latent_start_positive_analytical,atol=1e-12,rtol=4e-3))
    @test all(isapprox.(dMdt_latent_end_positive,dMdt_latent_end_positive_analytical,atol=1e-12,rtol=4e-3))
    @test all(isapprox.(dMdt_latent_start_negative,dMdt_latent_start_negative_analytical,atol=1e-12,rtol=4e-3))
    @test all(isapprox.(dMdt_latent_end_negative,dMdt_latent_end_negative_analytical,atol=1e-12,rtol=4e-3))

    dMdt_latent = dMdt_latent_start .+ dMdt_latent_end

    @test all(isapprox.(uu_test2[4*numofslugs+1:5*numofslugs],dMdt_latent,rtol=4e-3)) #dMdt

    F_start = ρₗ .* getδarea.(Ac,d,δstart)
    F_end = ρₗ .* getδarea.(Ac,d,δend)
    dLdt_start = -(ηplusvalue.*dMdt_latent_start_positive .+ ηminusvalue.*dMdt_latent_start_negative) ./ F_start
    dLdt_end = -(ηplusvalue.*dMdt_latent_end_positive .+ ηminusvalue.*dMdt_latent_end_negative) ./ F_end

    @test all(isapprox.(uu_test2[7*numofslugs+1:8*numofslugs],dLdt_start,atol=1e-12))
    @test all(isapprox.(uu_test2[8*numofslugs+1:9*numofslugs],dLdt_end,atol=1e-12))

    peri_start = peri .* (d .- 2δstart)/d
    peri_end = peri .* (d .- 2δend)/d

    C_start = ρₗ .* peri_start
    C_end = ρₗ .* peri_end
    dδdt_start = -(dMdt_latent_start .+ dLdt_start.*F_start) ./ C_start ./ Lfilm_start
    dδdt_end = -(dMdt_latent_end .+ dLdt_end.*F_end) ./ C_end ./ Lfilm_end

    @test all(isapprox.(uu_test2[5*numofslugs+1:6*numofslugs], dδdt_start,atol=1e-12))
    @test all(isapprox.(uu_test2[6*numofslugs+1:7*numofslugs], dδdt_end,atol=1e-12))

end

# test 3, five slugs, with heat transfer, given different initial film lengthes and liquid velocities to represent five different film states
@testset "five possible film states" begin

    ηplusvalue = rand()
    ηminusvalue = 0.0
    numofslugs = 5
    sys_tube = initialize_ohpsys(sys_plate,p_fluid,power,slugnum=numofslugs,ηplus=ηplusvalue)

    sys_tube.vapor.δstart .= zeros(numofslugs) .+ 1e-5 .+ 1e-5 .* rand(numofslugs) # initial velocity of the slugs
    sys_tube.vapor.δend .= zeros(numofslugs) .+ 1e-5 .+ 1e-5 .* rand(numofslugs) # initial velocity of the slugs
    sys_tube.liquid.dXdt[1:5] = [zero.(X) .- 1.0 .* rand() for X in sys_tube.liquid.dXdt[1:5]]# initial velocity of the slugs
   


    closedornot = sys_tube.tube.closedornot
    σ = sys_tube.liquid.σ
    L = sys_tube.tube.L
    ρₗ = sys_tube.liquid.ρ
    μₗ = sys_tube.liquid.μₗ
    ad_fac = sys_tube.vapor.ad_fac
    d = sys_tube.tube.d
    Ac = sys_tube.tube.Ac
    Xp = sys_tube.liquid.Xp
    Hfg = sys_tube.tube.PtoHfg.(sys_tube.vapor.P)
    k = sys_tube.vapor.k
    peri = sys_tube.tube.peri
    δstart = sys_tube.vapor.δstart
    δend = sys_tube.vapor.δend
    Lfilm_start = sys_tube.vapor.Lfilm_start
    Lfilm_end = sys_tube.vapor.Lfilm_end

    Lvaporplug = XptoLvaporplug(Xp,L,closedornot)

    # set up five different film states
    V = [mean(elem) for elem in sys_tube.liquid.dXdt]
    Vavg = mean(abs.(V))
    Ca = getCa.(μₗ,σ,Vavg)

    δdep = Catoδ(d,Ca,adjust_factor=ad_fac)
    Adeposit = getAdeposit(sys_tube,δdep)
    Adeposit_left = [elem[1] for elem in Adeposit]
    Adeposit_right = [elem[2] for elem in Adeposit]
    V_normal_start = circshift(Ac ./ (Ac .- Adeposit_right) .*  V,1)
    V_normal_end = Ac ./ (Ac .- Adeposit_left) .*  V

    Astart = getδarea(Ac,d,δstart)
    Aend = getδarea(Ac,d,δend)
    V_start_case5 = Ac / (Ac - Aend[5]) *  V[5]
    V_end_case5   = Ac / (Ac - Astart[4]) *  V[4]


    #case 2 to 5: set different initial film thicknesses
    Lfilm_start[2] = 0.1*sys_tube.wall.L_newbubble
    Lfilm_start[3] = Lvaporplug[3]/2 .- 0.1*sys_tube.wall.L_newbubble
    Lfilm_start[4] = Lvaporplug[4] .- 0.2*sys_tube.wall.L_newbubble
    Lfilm_start[5] = 0.1*sys_tube.wall.L_newbubble

    Lfilm_end[2] = 0.1*sys_tube.wall.L_newbubble
    Lfilm_end[3] = Lvaporplug[3]/2 .- 0.1*sys_tube.wall.L_newbubble
    Lfilm_end[4] = 0.1*sys_tube.wall.L_newbubble
    Lfilm_end[5] = Lvaporplug[5] .- 0.2*sys_tube.wall.L_newbubble

    u_tube = newstate(sys_tube) # initialize OHP tube

    uu_test3 = dynamicsmodel(u_tube[1:9*numofslugs],sys_tube)
    # dX2,ddXdt2,dM2,dδstart2,dδend2,dLfilm_start2,dLfilm_end2 = ComputationalHeatTransfer.vectoXMδL(uu_test2)

    V_start_analytical = [V_normal_start[1],
                          V_normal_start[2],
                          V_normal_start[3],
                          V_normal_start[4],
                          V_normal_start[5]]
    V_end_analytical = [V_normal_end[1],
                        V[2],
                        V_normal_end[3],
                        V_end_case5,
                        V_normal_end[5]]

    @test all(isapprox.(uu_test3[1:2:2*numofslugs-1],V_end_analytical,atol=1e-12))
    @test all(isapprox.(uu_test3[2:2:2*numofslugs],circshift(V_start_analytical,-1),atol=1e-12))


end

# test 4, 2 slugs, with heat transfer, given initial temperature differences, with initial velocity, test mass conservation δ L
@testset "mass conservation δ L" begin

    ηplusvalue = rand()
    ηminusvalue = 0.0
    numofslugs = 2
    sys_tube = initialize_ohpsys(sys_plate,p_fluid,power,slugnum=numofslugs,ηplus=ηplusvalue)
    sys_tube.liquid.dXdt = [zero.(X) .+ 1.0 .* rand() for X in sys_tube.liquid.dXdt]# initial velocity of the slugs
   
    ΔT = rand()
    ΔT_array = [ΔT,-ΔT]
    sys_tube.vapor.P[1] = sys_tube.tube.TtoP(Tref-ΔT_array[1])
    sys_tube.vapor.P[2] = sys_tube.tube.TtoP(Tref-ΔT_array[2])

    closedornot = sys_tube.tube.closedornot
    P1 = sys_tube.vapor.P[1]
    P2 = sys_tube.vapor.P[2]
    σ = sys_tube.liquid.σ
    L = sys_tube.tube.L
    ρₗ = sys_tube.liquid.ρ
    μₗ = sys_tube.liquid.μₗ
    ad_fac = sys_tube.vapor.ad_fac
    d = sys_tube.tube.d
    Ac = sys_tube.tube.Ac
    Xp = sys_tube.liquid.Xp
    Hfg = sys_tube.tube.PtoHfg.(sys_tube.vapor.P)
    k = sys_tube.vapor.k
    peri = sys_tube.tube.peri
    δstart = sys_tube.vapor.δstart
    δend = sys_tube.vapor.δend
    Lfilm_start = sys_tube.vapor.Lfilm_start
    Lfilm_end = sys_tube.vapor.Lfilm_end

    u_tube = newstate(sys_tube) # initialize OHP tube

    uu_test2 = dynamicsmodel(u_tube[1:9*numofslugs],sys_tube)


    V = [mean(elem) for elem in sys_tube.liquid.dXdt]
    # get a characteristic Capilarry number based on the average velocities
    Vavg = mean(abs.(V))
    Ca = getCa.(μₗ,σ,Vavg)

    δdep = Catoδ(d,Ca,adjust_factor=ad_fac)
    Adeposit = getAdeposit(sys_tube,δdep)
    Adeposit_left = [elem[1] for elem in Adeposit]
    Adeposit_right = [elem[2] for elem in Adeposit]
    Adeposit_start = circshift(Adeposit_right,1)
    Adeposit_end = Adeposit_left

    # run the dynamics model to get numerical dMdt_latent
    Xpvapor = getXpvapor(Xp,L,closedornot)
    dMdt_latent_start,dMdt_latent_end,dMdt_latent_start_positive,dMdt_latent_end_positive = dMdtdynamicsmodel(Xpvapor,sys_tube)
    dMdt_latent_start_negative = dMdt_latent_start .- dMdt_latent_start_positive
    dMdt_latent_end_negative = dMdt_latent_end .- dMdt_latent_end_positive

    dMdt_latent = dMdt_latent_start .+ dMdt_latent_end

    @test all(isapprox.(uu_test2[4*numofslugs+1:5*numofslugs],dMdt_latent,rtol=4e-3)) #dMdt

    F_start = ρₗ .* getδarea.(Ac,d,δstart)
    F_end = ρₗ .* getδarea.(Ac,d,δend)
    dLdt_start_analytical = -(ηplusvalue.*dMdt_latent_start_positive .+ ηminusvalue.*dMdt_latent_start_negative) ./ F_start .- Ac ./ (Ac .- Adeposit_start) .*  circshift(V,-1)
    dLdt_end_analytical  = -(ηplusvalue.*dMdt_latent_end_positive .+ ηminusvalue.*dMdt_latent_end_negative) ./ F_end .+ Ac ./ (Ac .- Adeposit_end) .*  V

    dLdt_start = uu_test2[7*numofslugs+1:8*numofslugs]
    dLdt_end = uu_test2[8*numofslugs+1:9*numofslugs]
    @test all(isapprox.(dLdt_start_analytical,dLdt_start,atol=1e-12))
    @test all(isapprox.(dLdt_end_analytical,dLdt_end,atol=1e-12))

    peri_start = peri .* (d .- 2δstart)/d
    peri_end = peri .* (d .- 2δend)/d

    C_start = ρₗ .* peri_start
    C_end = ρₗ .* peri_end
    dδdt_start = -(dMdt_latent_start .+ dLdt_start.*F_start) ./ C_start ./ Lfilm_start
    dδdt_end = -(dMdt_latent_end .+ dLdt_end.*F_end) ./ C_end ./ Lfilm_end

    dδdt_start = uu_test2[5*numofslugs+1:6*numofslugs]
    dδdt_end = uu_test2[6*numofslugs+1:7*numofslugs]
    # @test all(isapprox.(uu_test2[5*numofslugs+1:6*numofslugs], dδdt_start,atol=1e-12))
    # @test all(isapprox.(uu_test2[6*numofslugs+1:7*numofslugs], dδdt_end,atol=1e-12))

    dMdt_start_lhs = F_start .* dLdt_start .+ C_start .* Lfilm_start .* dδdt_start
    dMdt_end_lhs = F_end .* dLdt_end .+ C_end .* Lfilm_end .* dδdt_end
    dMdt_start_rhs = -dMdt_latent_start - ρₗ .* Adeposit_start .* Ac ./ (Ac .- Adeposit_start) .*  circshift(V,-1)
    dMdt_end_rhs = -dMdt_latent_end + ρₗ .* Adeposit_end .* Ac ./ (Ac .- Adeposit_end) .*  V

    @test all(isapprox.(dMdt_start_lhs,dMdt_start_rhs,atol=1e-12))
    @test all(isapprox.(dMdt_end_lhs,dMdt_end_rhs,atol=1e-12))
end

@testset "liquidmodel" begin

    # p_fluid = SaturationFluidProperty("Butane",Tref)

    ηplusvalue = rand()
    ηminusvalue = 0.0
    numofslugs = 2
    sys_tube = initialize_ohpsys(sys_plate,p_fluid,power,slugnum=numofslugs,ηplus=ηplusvalue)

    closedornot = sys_tube.tube.closedornot

    N_varyT = 10
    dTs = [1.0,-1.0]

    sys_tube.liquid.θarrays[1][N_varyT] += dTs[1]
    sys_tube.liquid.θarrays[2][N_varyT] += dTs[2]

    # use the liquidmodel function to get the du
    dus = ComputationalHeatTransfer.liquidmodel(sys_tube)


    α = sys_tube.liquid.α
    # k = sys_tube.liquid.kₗ
    Cpₗ = sys_tube.liquid.Cp
    Hₗ = sys_tube.liquid.Hₗ
    ρₗ = sys_tube.liquid.ρ
    peri = sys_tube.tube.peri
    Ac = sys_tube.tube.Ac

    @test all(isapprox.(α,p_fluid.kₗ/(Cpₗ*ρₗ),atol=1e-12))

    H_rhs = peri / (ρₗ*Cpₗ*Ac)

    # hand calculate the du
    θarrays = sys_tube.liquid.θarrays
    dx1 = mod(sys_tube.liquid.Xarrays[1][2] - sys_tube.liquid.Xarrays[1][1], sys_tube.tube.L)

    dθarray1 = zero(θarrays[1])
    dθarray1[N_varyT-1] = α .* dTs[1] ./ dx1 ./ dx1
    dθarray1[N_varyT] = α .* -2dTs[1] ./ dx1 ./ dx1 - Hₗ .* dTs[1] .* peri ./ (ρₗ*Cpₗ*Ac)
    dθarray1[N_varyT+1] = α .* dTs[1] ./ dx1 ./ dx1

    @test all(isapprox.(dθarray1,dus[1],atol=1e-12))


    dx2 = mod(sys_tube.liquid.Xarrays[2][2] - sys_tube.liquid.Xarrays[2][1], sys_tube.tube.L)

    dθarray2 = zero(θarrays[2])
    dθarray2[N_varyT-1] = α .* dTs[2] ./ dx2 ./ dx2
    dθarray2[N_varyT] = α .* -2dTs[2] ./ dx2 ./ dx2 - Hₗ .* dTs[2] .* peri ./ (ρₗ*Cpₗ*Ac)
    dθarray2[N_varyT+1] = α .* dTs[2] ./ dx2 ./ dx2

    @test all(isapprox.(dθarray2,dus[2],atol=1e-12))

    u_tube = newstate(sys_tube) # initialize OHP tube
    


    liquiddu = duliquidθtovec(dus)
    liquiddu_analytical = [0.0; dθarray1;
                           0.0; dθarray2]


    @test all(isapprox.(liquiddu,liquiddu_analytical,atol=1e-12))    
end


@testset "liquid Nu and Hₗ" begin

    Nuₗ=3.6

    sys_tube = initialize_ohpsys(sys_plate,p_fluid,power,Nu=Nuₗ)

    Hₗ_analytical = Nuₗ * p_fluid.kₗ / sys_tube.tube.d

    @test all(isapprox.(sys_tube.liquid.Hₗ,Hₗ_analytical,atol=1e-12))    
end

@testset "wall heat flux" begin

    numofslugs = 2
    sys_tube = initialize_ohpsys(sys_plate,p_fluid,power,slugnum=numofslugs)
    sys_tube.wall.θarray .+= [rand() for i in 1:length(sys_tube.wall.θarray)]
    u_tube = newstate(sys_tube) # initialize OHP tube
    sys_tube = getcurrentsys!(u_tube,sys_tube)

    qwall = sys_to_heatflux(sys_tube)

    xs = sys_tube.wall.Xarray

    @test all(isapprox.(sys_tube.wall.θarray,sys_tube.mapping.θ_interp_walltoliquid(xs),atol=1e-12))


    qwall_analytical = sys_tube.mapping.H_interp_liquidtowall(xs) .* (sys_tube.mapping.θ_interp_walltoliquid(xs) .- sys_tube.mapping.θ_interp_liquidtowall(xs)) .* sys_tube.tube.peri

    @test all(isapprox.(qwall,qwall_analytical,atol=1e-12))    
end

