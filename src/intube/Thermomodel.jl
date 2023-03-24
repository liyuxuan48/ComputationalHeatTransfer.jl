using Statistics

export dMdtdynamicsmodel,dynamicsmodel_steadyfilm,wallmodel,liquidmodel,dynamicsmodel,sys_to_heatflux,sys_to_Harray,integrator_to_heatflux,integrator_to_Harray


function dynamicsmodel(u::Array{Float64,1},p::PHPSystem)

    du = zero(deepcopy(u))

    sys = deepcopy(p)

    σ = sys.liquid.σ
    μₗ = sys.liquid.μₗ
    Xp = sys.liquid.Xp
    numofliquidslug = length(Xp)
    d = sys.tube.d
    peri = sys.tube.peri
    Ac = sys.tube.Ac
    angle = sys.tube.angle
    g = sys.tube.g
    P = sys.vapor.P
    Eratio_plus = sys.vapor.Eratio_plus
    Eratio_minus = sys.vapor.Eratio_minus
    δstart = sys.vapor.δstart
    δend = sys.vapor.δend
    Lfilm_start = sys.vapor.Lfilm_start
    Lfilm_end = sys.vapor.Lfilm_end

    V = [elem[2] for elem in sys.liquid.dXdt]
    Vavg = mean(abs.(V))
    Ca = getCa.(μₗ,σ,Vavg)
    
    ad_fac = sys.vapor.ad_fac
    δdeposit = Catoδ(d,Ca,adjust_factor=ad_fac)
    Adeposit = getAdeposit(sys,δdeposit)
    Adeposit_left = [elem[1] for elem in Adeposit]
    Adeposit_right = [elem[2] for elem in Adeposit]
    # Adeposit_vapor_left = [Adeposit_right[end];Adeposit_right[1:end-1]]
    # Adeposit_vapor_right = Adeposit_left
    Astart = getδarea(Ac,d,δstart)
    Aend = getδarea(Ac,d,δend)

    Lvaporplug = XptoLvaporplug(Xp,sys.tube.L,sys.tube.closedornot)
    Lliquidslug = XptoLliquidslug(Xp,sys.tube.L)

    heightg_interp = sys.mapping.heightg_interp
    Xp1 = [elem[1] for elem in sys.liquid.Xp]
    Xp2 = [elem[2] for elem in sys.liquid.Xp]
    heightg = map(tuple,heightg_interp(Xp1),heightg_interp(Xp2))

    # println(heightg)

    Xpvapor = getXpvapor(Xp,sys.tube.L,sys.tube.closedornot)

    ρₗ = sys.liquid.ρ
# get differential equation factors
    lhs = ρₗ*Ac .* Lliquidslug
    rhs_press = Ac ./ lhs

    Re_list = ρₗ .* abs.(V) .* d ./ μₗ
    f_coefficient = f_churchill.(Re_list .+ 1e-4)
    dXdt_to_stress = -0.125 .* f_coefficient .* ρₗ .* V .* abs.(V)
    rhs_dXdt = peri .* Lliquidslug .* dXdt_to_stress ./ lhs
    rhs_g = Ac*ρₗ ./ lhs

    if sys.tube.closedornot == false
        println("open loop not supported!")
         return "open loop not supported!"
        end

    if sys.tube.closedornot == true

        numofvaporbubble = numofliquidslug

        v_liquid_left=zeros(numofliquidslug)
        v_liquid_right=zeros(numofliquidslug)
        v_vapor_left_normal=zeros(numofvaporbubble)
        v_vapor_right_normal=zeros(numofvaporbubble)

        v_momentum = u[2*numofliquidslug+1:2:4*numofliquidslug]
        v_momentum_vapor_end = v_momentum
        v_momentum_vapor_start = [v_momentum[end];v_momentum[1:end-1]]
        v_liquid_left  = v_momentum .+ v_momentum .* Adeposit_left ./ (Ac .- Adeposit_left)
        v_liquid_right = v_momentum .+ v_momentum .* Adeposit_right ./(Ac .- Adeposit_right)
        
        rhs_dLdt = -v_momentum .* (v_liquid_right .- v_liquid_left) ./ Lliquidslug

            # vapor δ
            v_vapor_left_normal[2:end] = v_liquid_right[1:end-1]
            v_vapor_left_normal[1] = v_liquid_right[end]
            v_vapor_right_normal = v_liquid_left

            A_dδdt_right_vapor = Adeposit_left
            A_dδdt_right_liquid = Adeposit_right
            A_dδdt_left_vapor = zeros(size(A_dδdt_right_vapor))
            A_dδdt_left_vapor[2:end] = A_dδdt_right_liquid[1:end-1]
            A_dδdt_left_vapor[1] = A_dδdt_right_liquid[end]

            dMdt_latent_start,dMdt_sensible,dMdt_latent_end = dMdtdynamicsmodel(Xpvapor,sys)
            dMdt_latent_start_positive,dMdt_sensible_positive,dMdt_latent_end_positive = dMdtdynamicsmodel_positive(Xpvapor,sys)
            dMdt_latent_start_negative = dMdt_latent_start .- dMdt_latent_start_positive
            dMdt_latent_end_negative =   dMdt_latent_end   .- dMdt_latent_end_positive

            F_start = ρₗ .* Ac .* 4 .* δstart .* (d .- δstart) ./ (d^2)
            C_start = ρₗ .* Ac .* 4 .* (d .- 2δstart) ./ (d^2)

            F_end = ρₗ .* Ac .* 4 .* δend .* (d .- δend) ./ (d^2)
            C_end = ρₗ .* Ac .* 4 .* (d .- 2δend) ./ (d^2)

            L0threshold_film = 0.05*sys.wall.L_newbubble
            L0threshold_pure_vapor = 0.3*sys.wall.L_newbubble

            dLdt_start_normal = -(dMdt_latent_start_positive .* Eratio_plus .+ dMdt_latent_start_negative .* Eratio_minus) ./ F_start - v_vapor_left_normal
            dLdt_end_normal = -(dMdt_latent_end_positive .* Eratio_plus .+ dMdt_latent_end_negative .* Eratio_minus) ./ F_end + v_vapor_right_normal

            he_start_short = Bool.(heaviside.(-Lfilm_start .+ L0threshold_film))
            he_end_short = Bool.(heaviside.(-Lfilm_end .+ L0threshold_film))
            he_start_positive = Bool.(heaviside.(dLdt_start_normal))
            he_end_positive = Bool.(heaviside.(dLdt_end_normal))
            he_meet= Bool.(heaviside.(-Lvaporplug .+ Lfilm_start .+ Lfilm_end .+ L0threshold_pure_vapor))


            # zero dLdt case
            case2_flag_start = Bool.((1 .- he_meet) .* he_start_short .* (1 .- he_start_positive))
            # two ends meet and both nonzero case
            case3_flag_start = Bool.(he_meet .* he_start_short .* he_start_positive .+ he_meet .* (1 .- he_start_short) .* (1 .- he_end_short))
            # two ends meet and other side zero case
            case4_flag_start = Bool.(he_meet .* (1 .- he_start_short) .* he_end_short .* (1 .- he_end_positive))
            # two ends meet and this side zero case
            case5_flag_start = Bool.(he_meet .* he_start_short .* (1 .- he_start_positive))
            # normal case
            case1_flag_start = Bool.((1 .- case3_flag_start) .* (1 .- case2_flag_start) .* (1 .- case4_flag_start) .* (1 .- case5_flag_start))
            
            # zero dLdt case
            case2_flag_end = Bool.((1 .- he_meet) .* he_end_short .* (1 .- he_end_positive))
            # two ends meet and both nonzero case
            case3_flag_end = Bool.(he_meet .* he_end_short .* he_end_positive .+ he_meet .* (1 .- he_end_short) .* (1 .- he_start_short))
            # two ends meet and other side zero case
            case4_flag_end = Bool.(he_meet .* (1 .- he_end_short) .* he_start_short .* (1 .- he_start_positive))
            # two ends meet and this side zero case
            case5_flag_end = Bool.(he_meet .* he_end_short .* (1 .- he_end_positive))            
            # normal case
            case1_flag_end = Bool.((1 .- case3_flag_end) .* (1 .- case2_flag_end) .* (1 .- case4_flag_end) .* (1 .- case5_flag_end))

            he_matrix_start = hcat([case1_flag_start';case2_flag_start';case3_flag_start';case4_flag_start';case5_flag_start'])
            he_matrix_end   = hcat([case1_flag_end';case2_flag_end';case3_flag_end';case4_flag_end';case5_flag_end'])

            v_vapor_left_case5 = v_momentum_vapor_start .+ v_momentum_vapor_start .* Aend ./ (Ac .- Aend)
            v_vapor_right_case5 = v_momentum_vapor_end .+ v_momentum_vapor_end .* Astart ./ (Ac .- Astart)

            V_vapor_matrix_start = hcat([v_vapor_left_normal';v_momentum_vapor_start';v_vapor_left_normal'; v_vapor_left_normal';v_vapor_left_case5'])
            V_vapor_matrix_end   = hcat([v_vapor_right_normal';v_momentum_vapor_end';v_vapor_right_normal';v_vapor_right_normal';v_vapor_right_case5'])

            v_vapor_start_final = sum(he_matrix_start .* V_vapor_matrix_start,dims=1)
            v_vapor_end_final = sum(he_matrix_end .* V_vapor_matrix_end,dims=1)

            v_liquid_left_final = v_vapor_end_final
            v_liquid_right_final = [v_vapor_start_final[2:end];v_vapor_start_final[1]]

            dLdt_matrix_start = hcat([dLdt_start_normal';0 .* dLdt_start_normal';-v_vapor_left_normal';(v_vapor_right_case5 .- v_vapor_left_normal)';0 .* dLdt_start_normal'])
            dLdt_matrix_end   = hcat([dLdt_end_normal';0 .* dLdt_end_normal';v_vapor_right_normal';(v_vapor_right_normal .- v_vapor_left_case5)';0 .* dLdt_end_normal'])
        
            dLdt_start = sum(he_matrix_start .* dLdt_matrix_start,dims=1)
            dLdt_end = sum(he_matrix_end .* dLdt_matrix_end,dims=1)

            dδdt_start_normal = (-dMdt_latent_start .- F_start .* dLdt_start' .- ρₗ .* A_dδdt_left_vapor  .* v_vapor_left_normal) ./ (C_start .* Lfilm_start) 
            dδdt_end_normal = (-dMdt_latent_end     .- F_end   .* dLdt_end'   .+ ρₗ .* A_dδdt_right_vapor .* v_vapor_right_normal) ./ (C_end .* Lfilm_end)

            dδdt_start_case4 = (-dMdt_latent_start .- F_start .* dLdt_start' .- ρₗ .* A_dδdt_left_vapor  .* v_vapor_left_normal .+ ρₗ .* Astart  .* v_vapor_right_case5) ./ (C_start .* Lfilm_start) 
            dδdt_end_case4 = (-dMdt_latent_end     .- F_end   .* dLdt_end'   .+ ρₗ .* A_dδdt_right_vapor .* v_vapor_right_normal .- ρₗ .* Aend  .* v_vapor_left_case5) ./ (C_end .* Lfilm_end)


            dδdt_matrix_start = hcat([dδdt_start_normal';0 .* dδdt_start_normal';dδdt_start_normal';dδdt_start_case4';0 .* dδdt_start_normal'])
            dδdt_matrix_end   = hcat([dδdt_end_normal';0 .* dδdt_end_normal';dδdt_end_normal';dδdt_end_case4';0 .* dLdt_end_normal'])

            # println(size(dδdt_matrix_start))
        
            dδdt_start = sum(he_matrix_start .* dδdt_matrix_start,dims=1)
            dδdt_end = sum(he_matrix_end .* dδdt_matrix_end,dims=1)


        for i = 1:numofliquidslug
                du[2*i-1] = v_liquid_left_final[i]
                du[2*i] = v_liquid_right_final[i]

                if i != numofliquidslug
                    du[2*numofliquidslug + 2*i-1] = rhs_dXdt[i] + rhs_g[i]*(heightg[i][1]-heightg[i][end]) + rhs_press[i] * (P[i]-P[i+1]) + rhs_dLdt[i]
                else
                    du[2*numofliquidslug + 2*i-1] = rhs_dXdt[i] + rhs_g[i]*(heightg[i][1]-heightg[i][end]) + rhs_press[i] * (P[i]-P[1]) + rhs_dLdt[i]
                end                

                du[2*numofliquidslug + 2*i] = du[2*numofliquidslug + 2*i-1]

            end

            du[4*numofliquidslug+1:5*numofliquidslug] = dMdt_latent_start+dMdt_latent_end
            du[5*numofliquidslug+1:6*numofliquidslug] = dδdt_start # equals to 0 for now
            du[6*numofliquidslug+1:7*numofliquidslug] = dδdt_end # equals to 0 for now
            du[7*numofliquidslug+1:8*numofliquidslug] = dLdt_start # equals to 0 for now
            du[8*numofliquidslug+1:9*numofliquidslug] = dLdt_end # equals to 0 for now

            return du
    end


end

function dMdtdynamicsmodel(Xpvapor::Array{Tuple{Float64,Float64},1},sys::PHPSystem)

    dMdt_sensible=zeros(length(Xpvapor))
    dMdt_latent_start=zeros(length(Xpvapor))
    dMdt_latent_end=zeros(length(Xpvapor))

    L = sys.tube.L

    peri = sys.tube.peri
    Ac = sys.tube.Ac
    k = sys.vapor.k

    Lfilm_start = sys.vapor.Lfilm_start
    Lfilm_end = sys.vapor.Lfilm_end

    P = sys.vapor.P
    θarrays = sys.liquid.θarrays
    Xarrays = sys.liquid.Xarrays

    @unpack PtoT,PtoHfg = sys.tube
    θ = PtoT.(P)

    Hfg = PtoHfg.(P)

    
    H_interp = sys.mapping.H_interp_liquidtowall
    θ_wall_interp = sys.mapping.θ_interp_walltoliquid

    dx_wall = sys.wall.Xarray[2]-sys.wall.Xarray[1]

    Lvapor = XptoLvaporplug(sys.liquid.Xp,sys.tube.L,sys.tube.closedornot)
    Lvapor_pure = max.(Lvapor - Lfilm_start - Lfilm_end,0.0)

    for i = 1:length(Xpvapor)
        Nstart = Int64(max(2 , div(Lfilm_start[i],dx_wall)))
        heatflux_start = quad_trap(H_interp,θ_wall_interp,θ[i],Xpvapor[i][1],mod(Xpvapor[i][1]+Lfilm_start[i],L),L,Nstart)

        Nvapor_pure = Int64(max(2 , div(Lvapor_pure[i],dx_wall)))
        heatflux_pure_vapor = quad_trap(H_interp,θ_wall_interp,θ[i],mod(Xpvapor[i][1]+Lfilm_start[i],L),mod(Xpvapor[i][2]-Lfilm_end[i],L),L,Nvapor_pure)
        Nend = Int64(max(2 , div(Lfilm_end[i],dx_wall)))
        heatflux_end = quad_trap(H_interp,θ_wall_interp,θ[i],mod(Xpvapor[i][2]-Lfilm_end[i],L),Xpvapor[i][2],L,Nend)

        slope_r = getslope(θarrays[i][2],θarrays[i][1],Xarrays[i][2],Xarrays[i][1])
        slope_l = (i == 1) ? getslope(θarrays[1][end],θarrays[1][end-1],Xarrays[1][end],Xarrays[1][end-1]) : getslope(θarrays[i-1][end],θarrays[i-1][end-1],Xarrays[i-1][end],Xarrays[i-1][end-1])


        axial_rhs_end = Ac*k*slope_r /Hfg[i]
        axial_rhs_start = Ac*k*(-slope_l) /Hfg[i]

        dMdt_latent_start[i] = heatflux_start*peri/Hfg[i] + axial_rhs_start

        dMdt_sensible[i] = heatflux_pure_vapor*peri/Hfg[i] .* 0.0

        dMdt_latent_end[i] = heatflux_end*peri/Hfg[i] + axial_rhs_end
        
        end
            return dMdt_latent_start,dMdt_sensible,dMdt_latent_end
end

function dMdtdynamicsmodel_positive(Xpvapor::Array{Tuple{Float64,Float64},1},sys::PHPSystem)

    dMdt_sensible=zeros(length(Xpvapor))
    dMdt_latent_start=zeros(length(Xpvapor))
    dMdt_latent_end=zeros(length(Xpvapor))

    L = sys.tube.L

    peri = sys.tube.peri
    Ac = sys.tube.Ac
    k = sys.vapor.k

    Lfilm_start = sys.vapor.Lfilm_start
    Lfilm_end = sys.vapor.Lfilm_end

    P = sys.vapor.P
    θarrays = sys.liquid.θarrays
    Xarrays = sys.liquid.Xarrays

    @unpack PtoT,PtoHfg = sys.tube
    θ = PtoT.(P)

    Hfg = PtoHfg.(P)

    
    H_interp = sys.mapping.H_interp_liquidtowall
    θ_wall_interp = sys.mapping.θ_interp_walltoliquid

    dx_wall = sys.wall.Xarray[2]-sys.wall.Xarray[1]

    Lvapor = XptoLvaporplug(sys.liquid.Xp,sys.tube.L,sys.tube.closedornot)
    Lvapor_pure = max.(Lvapor - Lfilm_start - Lfilm_end,0.0)

    for i = 1:length(Xpvapor)
        Nstart = Int64(max(2 , div(Lfilm_start[i],dx_wall)))
        heatflux_start = quad_trap_positive(H_interp,θ_wall_interp,θ[i],Xpvapor[i][1],mod(Xpvapor[i][1]+Lfilm_start[i],L),L,Nstart)

        Nvapor_pure = Int64(max(2 , div(Lvapor_pure[i],dx_wall)))
        heatflux_pure_vapor = quad_trap_positive(H_interp,θ_wall_interp,θ[i],mod(Xpvapor[i][1]+Lfilm_start[i],L),mod(Xpvapor[i][2]-Lfilm_end[i],L),L,Nvapor_pure)

        Nend = Int64(max(2 , div(Lfilm_end[i],dx_wall)))
        heatflux_end = quad_trap_positive(H_interp,θ_wall_interp,θ[i],mod(Xpvapor[i][2]-Lfilm_end[i],L),Xpvapor[i][2],L,Nend)

        slope_r = getslope(θarrays[i][2],θarrays[i][1],Xarrays[i][2],Xarrays[i][1])
        slope_l = (i == 1) ? getslope(θarrays[1][end],θarrays[1][end-1],Xarrays[1][end],Xarrays[1][end-1]) : getslope(θarrays[i-1][end],θarrays[i-1][end-1],Xarrays[i-1][end],Xarrays[i-1][end-1])


        axial_rhs_end = Ac*k*slope_r /Hfg[i]
        axial_rhs_start = Ac*k*(-slope_l) /Hfg[i]

        dMdt_latent_start[i] = heatflux_start*peri/Hfg[i] + axial_rhs_start

        dMdt_sensible[i] = heatflux_pure_vapor*peri/Hfg[i] .* 0.0

        dMdt_latent_end[i] = heatflux_end*peri/Hfg[i] + axial_rhs_end
        
        end
            return dMdt_latent_start,dMdt_sensible,dMdt_latent_end
end

function liquidmodel(p::PHPSystem)
    sys = deepcopy(p)
    θarrays = sys.liquid.θarrays
    # nondihv_tonondihl = 0.0046206704347650325 # temperary variable to fix different nondimensionlaization

    du = zero.(deepcopy(θarrays))

    # γ = sys.vapor.γ
    Ac = sys.tube.Ac

    Hₗ = sys.liquid.Hₗ
    peri = sys.tube.peri
    α = sys.liquid.α
    Cpₗ = sys.liquid.Cp
    ρₗ = sys.liquid.ρ

    H_rhs = peri / (ρₗ*Cpₗ*Ac)

    for i in eachindex(θarrays)
        
        xs = sys.liquid.Xarrays[i];
        dx = mod(xs[2] - xs[1], sys.tube.L)

        θ_wall_inter = sys.mapping.θ_interp_walltoliquid

        fx = map(θ_wall_inter, xs) - θarrays[i]
        du[i] = α .* laplacian(θarrays[i]) ./ dx ./ dx + Hₗ .* fx .* H_rhs
    end
    return du
end


"""
    This is a function to get the laplacian of a vector field u

    For now zero-gradient boundary condition is used.

    u    ::  an array
"""


function laplacian(u)
    unew = deepcopy(u)

    dl = ones(length(u)-1)
    dr = dl
    d  = -2*ones(length(u))
    d[1] = -1
    d[end]= -1

    A = Tridiagonal(dl, d, dr)

    unew = A*u

    return (unew)
end


# q'
function sys_to_heatflux(p::PHPSystem)

    sys = deepcopy(p)

    θarray = sys.wall.θarray
    # γ = sys.vapor.γ
    Hₗ = sys.liquid.Hₗ
    # He = sys.evaporator.He
    # k = sys.vapor.k
    # δ = sys.vapor.δ
    # Hvapor = k ./ δ

    peri = sys.tube.peri

    # dx = sys.wall.Xarray[2]-sys.wall.Xarray[1]

    # Xarray = sys.wall.Xarray
    θ_interp_liquidtowall = sys.mapping.θ_interp_liquidtowall
    H_interp_liquidtowall = sys.mapping.H_interp_liquidtowall

    xs =  sys.wall.Xarray

    dθarray = map(θ_interp_liquidtowall, xs) .- θarray
    Harray  = map(H_interp_liquidtowall, xs)

    # qwallarray = -Harray.*dθarray
    qwallarray = -Harray.*dθarray*peri

    # println(Harray)
end


function integrator_to_heatflux(inte)
    sys_to_heatflux(inte.p)
end

function sys_to_Harray(p::PHPSystem)

    sys = deepcopy(p)

    H_interp_liquidtowall = sys.mapping.H_interp_liquidtowall

    xs =  sys.wall.Xarray

    Harray  = map(H_interp_liquidtowall, xs)

    Harray
end

function quad_trap(H_interp,θ_interp,θvapor_one, a,b,L,N) 
    h = mod((b-a),L)/N
    
    # println(h)

    H_interp(a)*(θ_interp(a)-θvapor_one)
    # println(b)

    int = h * ( H_interp(a)*(θ_interp(a)-θvapor_one) + H_interp(b)*(θ_interp(b)-θvapor_one) ) / 2
    for k=1:N-1
        xk = mod((b-a),L) * k/N + a
        int = int + h*H_interp(mod(xk,L))*(θ_interp(mod(xk,L))-θvapor_one)
    end
    return (b>a) ? int : 0
end

function quad_trap_positive(H_interp,θ_interp,θvapor_one, a,b,L,N) 
    h = mod((b-a),L)/N
    
    # println(h)

    H_interp(a)*(θ_interp(a)-θvapor_one)
    # println(b)

    int = maximum([h * ( H_interp(a)*(θ_interp(a)-θvapor_one) + H_interp(b)*(θ_interp(b)-θvapor_one) ) / 2, 0.0])
    for k=1:N-1
        xk = mod((b-a),L) * k/N + a
        int = int + maximum([h*H_interp(mod(xk,L))*(θ_interp(mod(xk,L))-θvapor_one),0.0])
    end
    return (b>a) ? int : 0
end

function integrator_to_Harray(inte)
    sys_to_Harray(inte.p)
end

getslope(y2,y1,x2,x1) = (y2-y1)/(x2-x1)

heaviside(x::AbstractFloat) = ifelse(x < 0, zero(x), ifelse(x > 0, one(x), oftype(x,0.0)))