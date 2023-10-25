import ConstrainedSystems: init, solve

export ODE_innertube,ODE_steadyfilm,timemarching!

function ODE_innertube(u,p,t)

    sys_init = p

    index_dynamics_end = findfirst(x->abs(x+1e10) <= 10^(-1), u)

    newsys = getcurrentsys!(u,sys_init)

    # println(sys_init.tube.PtoT(sys_init.vapor.P))
    # println(newsys.tube.PtoT(newsys.vapor.P))

    dynamicsdu = dynamicsmodel(u[1:index_dynamics_end-1],newsys)

    liquiddu = duliquidθtovec(liquidmodel(newsys))

    du = [dynamicsdu;liquiddu]
    
    return(du)

end

function ODE_steadyfilm(u,p,t)

    index_dynamics_end = findfirst(x->abs(x+1e10) <= 10^(-1), u)

    newsys = getcurrentsys!(u,p)

    dynamicsdu = dynamicsmodel_steadyfilm(u[1:index_dynamics_end-1],newsys)

    liquiddu = duliquidθtovec(liquidmodel(newsys))

    du = [dynamicsdu;liquiddu]

    return(du)

end


# weakly coupled alternate time marching
function timemarching!(integrator_tube,integrator_plate,tstep::Float64)

    currentsys = getcurrentsys!(integrator_tube.u,integrator_tube.p)
    currentsys.wall.θarray = temperature_linesource(integrator_plate)

    qtmp = sys_to_heatflux(currentsys)
    sys_plate = integrator_plate.p
    set_linesource_strength!(sys_plate,qtmp)
    
    step!(integrator_tube,tstep,true);
    ADI_timemarching!(temperature(integrator_plate),sys_plate,tstep)
    integrator_plate.t += tstep
    integrator_tube,integrator_plate
end

function init(u_tube::Vector{Float64},tspan::Tuple{Any, Any},sys_tube::PHPSystem,kwargs...)
    
    cb_boiling =  DiscreteCallback(boiling_condition,boiling_affect!)
    cb_vapormerging =  DiscreteCallback(merging_condition,merging_affect!)
    cb_liquidmerging = DiscreteCallback(vaporMergingCondition,vaporMergingAffect!)
    cb_fixdx =  DiscreteCallback(fixdx_condition,fixdx_affect!)
    cb_slugbc = DiscreteCallback(slugbc_condition,slugbc_affect!)
    cbst = CallbackSet(cb_fixdx,cb_boiling,cb_vapormerging,cb_liquidmerging,cb_slugbc);

    prob = ODEProblem(ODE_innertube, u_tube, tspan, sys_tube,kwargs...) # construct integrator_tube problem
    return init(prob, alg=RK4(),dt=1e-3,save_on=false, callback=cbst,maxiters=1e10,kwargs...)
end