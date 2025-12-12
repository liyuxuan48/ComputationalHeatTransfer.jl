using Statistics
using ComputationalHeatTransfer

# test Ca to δ correlation (Aussillous and Quere, 2000)
@testset "liquid film thickness correlation" begin
    Ns = 30
    d = 1e-3
    δmin = 2e-6
    δmax = 1e-4
    Ca = rand(Ns) .* 0.05
    δ = 1.34 * (Ca .^ (2/3)) ./ (1 .+ 3.35 .* (Ca .^ (2/3))) .* d/2
    δ_clamp = clamp.(δ,δmin,δmax)

    @test all(δ_clamp .≈ Catoδ.([d],Ca,δmin=δmin,δmax=δmax))

    ad_fac = 1 + 0.3 * rand()

    δ_clamp_adfac = clamp.(δ .* ad_fac,δmin,δmax)
    @test all(δ_clamp_adfac .≈ Catoδ.([d],Ca,adjust_factor=ad_fac,δmin=δmin,δmax=δmax))
end

# test Churchill friction factor correlation (Churchill, 1977) 
@testset "churchill friction factor" begin
    Ns = 30
    ϵ = 0.001 # relative roughness
    Re_laminar = rand(Ns)*200
    f_ch = f_churchill.(Re_laminar,[ϵ])

    @test all(isapprox.(f_ch, 64 ./ Re_laminar,atol=1e-12,rtol=1e-4))

    Re_turbulent = rand(Ns) .* 1e6 .+ 4000
    f_ch = f_churchill.(Re_turbulent,[ϵ])

    A = (-2.457 .* log.((7 ./ Re_turbulent) .^0.9 .+ 0.27 .* ϵ)) .^16
    B = (37530 ./ Re_turbulent) .^16
    f_ch_analytical = 8*((8 ./ Re_turbulent) .^12 .+ (1 ./ ((A .+ B) .^1.5))) .^ (1/12)

    @test all(isapprox.(f_ch,f_ch_analytical,atol=1e-12,rtol=1e-12))

end

# test boiling ΔT correlation (Qu and Ma, 2007)
@testset "boiling ΔT" begin
    Rn = 3e-6 + rand()*1e-6
    Tref = 273.15 + rand()*50
    d = 1e-3 + rand()*1e-3
    fluid_type = "Butane"
    PtoT,TtoP,PtoD,DtoP,PtoHfg = createCoolPropinterpolation(fluid_type)

    ΔT = RntoΔT(Rn,Tref,fluid_type,d,TtoP)


    p_fluid = SaturationFluidProperty(fluid_type,Tref)
    Rkg = p_fluid.R/p_fluid.M
    Rin = d/2
    P = TtoP(Tref)
    Hfg = p_fluid.hᵥ - p_fluid.hₗ

    y = Rkg .* Tref ./ Hfg .* log.(1 .+ 2 .* p_fluid.σ ./ P .* (1 ./ Rn - 1 ./ (2 .* (Rin - 0.0))))
    ΔT_analytical = Tref .* (1 ./ (1 .- y) - 1)

    @test isapprox.(ΔT,ΔT_analytical,atol=1e-12,rtol=1e-12)
end