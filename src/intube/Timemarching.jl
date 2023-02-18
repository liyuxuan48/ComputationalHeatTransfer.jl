export ODE_innertube,ODE_steadyfilm,timemarching!

function ODE_innertube(u,p,t)

    sys_init = p

    index_dynamics_end = findfirst(x->abs(x+1e10) <= 10^(-1), u)

    # println(u[1:100])

    newsys = getcurrentsys(u,sys_init)

    dynamicsdu = dynamicsmodel(u[1:index_dynamics_end-1],newsys)

    liquiddu = duliquidθtovec(liquidmodel(newsys))

    du = [dynamicsdu;liquiddu]

    return(du)

end

function ODE_steadyfilm(u,p,t)

    index_dynamics_end = findfirst(x->abs(x+1e10) <= 10^(-1), u)

    newsys = getcurrentsys(u,p)

    dynamicsdu = dynamicsmodel_steadyfilm(u[1:index_dynamics_end-1],newsys)

    liquiddu = duliquidθtovec(liquidmodel(newsys))

    du = [dynamicsdu;liquiddu]

    return(du)

end


# weakly coupled alternate time marching
function timemarching!(integrator_tube,integrator_plate,tstep::Float64)
    
    step!(integrator_tube,tstep,true);

    currentsys = deepcopy(getcurrentsys(integrator_tube.u,integrator_tube.p))
    currentsys.wall.θarray = deepcopy(temperature_linesource(integrator_plate))
    integrator_tube.p = deepcopy(currentsys)
    qtmp = deepcopy(sys_to_heatflux(currentsys))
    sys_plate = integrator_plate.p
    set_linesource_strength!(sys_plate,qtmp)

    ADI_timemarching!(temperature(integrator_plate),sys_plate,tstep)
    integrator_plate.t += tstep
    
    integrator_tube,integrator_plate
end