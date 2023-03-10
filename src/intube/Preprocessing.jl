export SimulationResult,randomXp,L_to_boiltime,initialize_ohpsys,store!

function randomXp(tube::Tube,numofslugs=30,chargeratio=0.46,σ_charge=0.1)

    L = tube.L
    Lmin = tube.d

    σ_persection = σ_charge*L/sqrt(numofslugs)

    L_perslug=L/numofslugs*chargeratio
    L_persection=L/numofslugs

    Ls = abs.((rand(numofslugs) .- 0.5).*σ_persection .+ L_perslug)

    Xp1s = zeros(numofslugs);
    Xp2s = deepcopy(Xp1s);

    if minimum(Ls) > Lmin && maximum(Ls) < L_persection

        for i in eachindex(Xp1s)
            Xp1s[i] = (i-1)*L_persection
            Xp2s[i] = Xp1s[i] + Ls[i]
        end

        displacement = L*rand()

        Xp1s = mod.(Xp1s.+displacement,L)
        Xp2s = mod.(Xp2s.+displacement,L)

        else println("generation failed")
    end

    X0 = map(tuple,Xp1s,Xp2s)
    real_ratio = sum(Ls)/L

    X0,real_ratio
end

function L_to_boiltime(L_newbubble,Rn,fluid_type,vapor::Vapor,tube::Tube)
    P = vapor.P[1]
    property = SaturationFluidProperty(fluid_type,PtoT(P));

    tube_d = tube.d
    ρₗ = property.ρₗ
    Cpₗ = property.Cpₗ 
    kₗ = property.kₗ

    Hfg = PtoHfg(P)
    Tref = PtoT(P)
    ρv = PtoD(P)
    ΔTthres = RntoΔT(Rn,Tref,fluid_type,tube_d)

    A = (2 * (ΔTthres/Tref) * Hfg*ρv/ρₗ)^0.5
    Ja = ΔTthres*Cpₗ*ρₗ/ρv/Hfg
    B = (12*kₗ/pi/ρₗ/Cpₗ)^0.5 * Ja

    t = 1e-3:1e-3:10e0
    tstar = t .* A^2 ./ B^2 

    Rplus = 2/3 .* ((tstar .+ 1).^1.5 .- (tstar).^1.5 .- 1);
    R = Rplus .* B^2 ./ A; 
    L_equivalent = (4/3) ./ tube_d .* R.^2 

    interp_Rtot = LinearInterpolation(L_equivalent, t);

    return interp_Rtot(L_newbubble)
end

function initialize_ohpsys(OHPtype,fluid_type,sys,p_fluid,Tref,power)

    L = (sys.qline[1].arccoord[1] + sys.qline[1].arccoord[end])  # total length of the pipe when streched to a 1D pipe (an approximate here)
    ohp = sys.qline[1].body
    N=numpts(ohp)

    # if OHPtype == "ASETS-II OHP 1 LARGE HEATER" || "ASETS-II OHP 2 LARGE HEATER" || "ASETS-II OHP 3 LARGE HEATER" || "ASETS-II OHP 1 SMALL HEATER" || "ASETS-II OHP 2 SMALL HEATER" || "ASETS-II OHP 3 SMALL HEATER" 
        # tube geometries
        tube_d = 1e-3; # tube diameter
        peri = 4*tube_d # Tube perimeter
        Ac = tube_d*tube_d # tube cross-sectional area
        L2D = 133.83*1e-3 # the actual length of the bended pipe in the real world
        angle = 0*pi/2 # inclination angle 
        closedornot = true

        tube = Tube(tube_d,peri,Ac,L,L2D,angle,gravity,closedornot,N,fluid_type);

        # liquid
        Nu = 3.60
        Hₗ = p_fluid.kₗ/tube_d * Nu # Nusselt number 3.60

        X0,realratio = randomXp(tube,30)
        dXdt0_l = zeros(length(X0))
        dXdt0_r = zeros(length(X0))
        dXdt0 = map(tuple,dXdt0_l,dXdt0_r);

        N=numpts(ohp)
        Xarrays,θarrays = constructXarrays(X0,N,Tref,L);

        liquids=Liquid(Hₗ,p_fluid.ρₗ,p_fluid.Cpₗ,p_fluid.αₗ,p_fluid.μₗ,p_fluid.σ,X0,dXdt0,Xarrays,θarrays);

        # Vapor
        Hᵥ = 0.0 # Nusselt number 4.36
        P_initial = 0*zeros(length(X0)) .+ TtoP(Tref);
        δfilm = 2e-5
        δstart_initial = 0*zeros(length(X0)) .+ δfilm ;
        δend_initial = 0*zeros(length(X0)) .+ δfilm ;

        Lvaporplug = XptoLvaporplug(X0,L,tube.closedornot)
        Lfilm_start_initial = 0.01 .* Lvaporplug
        Lfilm_end_initial = 0.01 .* Lvaporplug
        δmin = 2e-6
        # vapors=Vapor(ad_fac,Hᵥ,p_fluid.kₗ,δmin,Eratio_plus,Eratio_minus,P,δfilm_deposit,δstart,δend,Lfilm_start,Lfilm_end);
        vapors=Vapor(k = p_fluid.kₗ,P=P_initial,δstart=δstart_initial,δend=δend_initial,Lfilm_start=Lfilm_start_initial,Lfilm_end=Lfilm_end_initial);

        # Wall
        # ΔTthres = RntoΔT(Rn,Tref,fluid_type,tube_d)
        nucleatenum = 1000
        Xstations = sort(rand(nucleatenum) .* L);
        Xstation_time = zeros(nucleatenum);

        boil_type = "wall T"
        L_newbubble = 8tube_d
        # boil_interval = L_to_boiltime(L_newbubble,Rn,fluid_type,vapors::Vapor,tube::Tube)
        boil_interval = 1.0
        Xwallarray,θwallarray = constructXarrays(sys.qline[1].arccoord,L,Tref);
        θwallarray .= Tref

        wall = Wall(fluid_type=fluid_type,boil_type=boil_type,power=power,L_newbubble=L_newbubble,Xstations=Xstations,boiltime_stations=Xstation_time,Xarray=Xwallarray,θarray=θwallarray);
        # wall = Wall(fluid_type,boil_type,boil_interval,Rn,L_newbubble,Xstations,Xstation_time,Xwallarray,θwallarray);

    # end

    sys0_nomapping = PHPSystem_nomapping(tube,liquids,vapors,wall);
    θ_interp_walltoliquid, θ_interp_liquidtowall, H_interp_liquidtowall, P_interp_liquidtowall = sys_interpolation(sys0_nomapping)
    mapping = Mapping(θ_interp_walltoliquid, θ_interp_liquidtowall, H_interp_liquidtowall, P_interp_liquidtowall);

    sys0 = PHPSystem(tube,liquids,vapors,wall,mapping);

    # Lvaporplug = XptoLvaporplug(X0,sys0.tube.L,sys0.tube.closedornot)
    # M = PtoD(P) .* Lvaporplug .* Ac

    # u=[XMδLtovec(X0,dXdt0,M,δstart,δend,Lfilm_start,Lfilm_end); liquidθtovec(sys0.liquid.θarrays)];
    sys0
end

function newstate(sys0::PHPSystem)
    Ac = sys0.tube.Ac

    X0 = sys0.liquid.Xp
    dXdt0 = sys0.liquid.dXdt

    δstart = sys0.vapor.δstart
    δend = sys0.vapor.δend
    Lfilm_start = sys0.vapor.Lfilm_start
    Lfilm_end = sys0.vapor.Lfilm_end
    P = sys0.vapor.P
    

    Lvaporplug = XptoLvaporplug(X0,sys0.tube.L,sys0.tube.closedornot)
    M = PtoD(P) .* Lvaporplug .* Ac

    u=[XMδLtovec(X0,dXdt0,M,δstart,δend,Lfilm_start,Lfilm_end); liquidθtovec(sys0.liquid.θarrays)];
end

mutable struct SimulationResult
    tube_hist_t      ::Vector{Any}
    tube_hist_u      ::Vector{Any}
    tube_hist_θwall  ::Vector{Any}
    boil_hist        ::Vector{Any}
    plate_T_hist     ::Vector{Any}
    integrator_tube  ::Any
    integrator_plate ::Any
    grid             ::Any
end

function SimulationResult(int_tube,int_plate)
    
    boil_hist=[]
    plate_T_hist = []
    tube_hist_u  = []
    tube_hist_t = []
    tube_hist_θwall = []
    integrator_tube = deepcopy(int_tube)
    integrator_plate = deepcopy(int_plate)
    grid = int_plate.p.grid
    
    return SimulationResult(tube_hist_t,tube_hist_u,tube_hist_θwall,boil_hist,plate_T_hist,integrator_tube,integrator_plate,grid)
end

function store!(sr,integrator_tube,integrator_plate)
        
#         sr.boil_hist = deepcopy(boil_hist);
        push!(sr.plate_T_hist,deepcopy(temperature(integrator_plate)));
        push!(sr.tube_hist_θwall,deepcopy(integrator_tube.p.wall.θarray))
        push!(sr.tube_hist_u,deepcopy(integrator_tube.u));
        push!(sr.tube_hist_t,deepcopy(integrator_tube.t));
end