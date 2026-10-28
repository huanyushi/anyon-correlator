using JLD2, Printf

function lbfgs_direction(g,pairs)
    valid = filter(p->dot(p[1],p[2])>1e-10*norm(p[1])*norm(p[2]),pairs)
    isempty(valid) && return -g
    q = copy(g)
    alpha = zeros(length(valid))
    rho = [1/dot(s,y) for (s,y) in valid]
    for i in reverse(eachindex(valid))
        s,y = valid[i]
        alpha[i] = rho[i]*dot(s,q)
        q -= alpha[i]*y
    end
    s,y = last(valid)
    r = (dot(s,y)/dot(y,y))*q
    for i in eachindex(valid)
        s,y = valid[i]
        r += s*(alpha[i]-rho[i]*dot(y,r))
    end
    return dot(g,-r)<0 ? -r : -g
end

function write_step(outdir,row)
    open(joinpath(outdir,"trajectory.jsonl"),"a") do io
        write(io,JSON.json(row));write(io,'\n');flush(io)
    end
end

function optimize(x,D,hx,hz;chi=8,chi_check=16,maxsteps=200,budget=300,tol=1e-9,outdir)
    H,ops = operators(hx,hz)
    pairs = Tuple{Vector{Float64},Vector{Float64}}[]
    start = time()
    state = contract(x,D,H;chi,tol)
    measurement = measure(state.psi,state.env,ops,hx,hz)
    physical_ok(measurement,hx,hz) || error("initial contraction violates physical bounds")
    jldsave(joinpath(outdir,"initial.jld2");parameters=x,D,hx,hz,chi)
    gprev = nothing
    sprev = nothing
    status = "budget"
    gradnorm = Inf
    steps = 0
    for step in 0:maxsteps
        energy,g = energy_gradient(x,D,H,state;tol)
        gradnorm = norm(g)
        if !isnothing(gprev)
            y = g-gprev
            if dot(sprev,y)>1e-10*norm(sprev)*norm(y)
                push!(pairs,(sprev,y));length(pairs)>8 && popfirst!(pairs)
            end
        end
        measurement = measure(state.psi,state.env,ops,hx,hz)
        abs(measurement.energy-energy)<=1e-8 || error("Hamiltonian and observable energies disagree")
        row = (;step,energy,gradnorm,chi,residual=Float64(state.info.convergence_error),seconds=time()-start,measurement)
        write_step(outdir,row)
        @printf("AD step=%d e=%+.12f |g|=%.3e mx=%.7f mz=%.7f residual=%.2e wall=%.1fs\n",step,energy,gradnorm,measurement.mx,measurement.mz,state.info.convergence_error,time()-start);flush(stdout)
        jldsave(joinpath(outdir,"latest.jld2");parameters=x,D,hx,hz,chi,step,energy,gradnorm)
        if gradnorm<=1e-5
            status="stationary";break
        end
        (step==maxsteps || time()-start>=budget) && break
        direction = lbfgs_direction(g,pairs)
        alpha = min(1.0,0.2/norm(direction))
        accepted = false
        for attempt in 1:16
            trialx = x+alpha*direction
            evaluation = try
                trial = contract(trialx,D,H;chi,tol)
                tm = measure(trial.psi,trial.env,ops,hx,hz)
                objective_ok = physical_ok(tm,hx,hz) && trial.energy<=energy+1e-4*alpha*dot(g,direction)
                if objective_ok
                    check = contract(trialx,D,H;chi,tol,seed=1)
                    cm = measure(check.psi,check.env,ops,hx,hz)
                    objective_ok = physical_ok(cm,hx,hz) && abs(check.energy-trial.energy)<=1e-6 && check.energy<=energy+1e-7
                    if objective_ok
                        observable_delta = maximum(abs(tm.observables[k]-cm.observables[k]) for k in keys(tm.observables))
                        objective_ok &= observable_delta<=1e-5
                    end
                end
                (;ok=objective_ok,state=trial)
            catch err
                err isa ContractionFailure || rethrow()
                write_step(outdir,(;step,attempt,alpha,event="trial_rejected",error=sprint(showerror,err)))
                (;ok=false,state=nothing)
            end
            if evaluation.ok
                gprev,sprev = g,trialx-x
                x,state = trialx,evaluation.state
                steps=step+1
                jldsave(joinpath(outdir,@sprintf("step%04d.jld2",steps));parameters=x,D,hx,hz,chi,step=steps,alpha,energy=state.energy)
                accepted=true;break
            end
            alpha/=2
            time()-start>=budget && break
        end
        if !accepted
            status = time()-start>=budget ? "budget" : "armijo_failed"
            break
        end
    end
    # A final gradient on the actual frozen tensor; accepted-step g is pre-step.
    _,gfinal = energy_gradient(x,D,H,state;tol)
    gradnorm = norm(gfinal)
    status=="budget" && gradnorm<=1e-5 && (status="stationary")
    final_measurement = measure(state.psi,state.env,ops,hx,hz)
    jldsave(joinpath(outdir,"final.jld2");parameters=x,D,hx,hz,chi,steps,status,gradnorm,measurement=final_measurement)
    audit_result = audit(x,D,hx,hz;chi,chi_check,tol,gradnorm,status)
    summary = (;status,steps,gradnorm,seconds=time()-start,measurement=final_measurement,reference=reference_values(hx,hz),audit=audit_result)
    save_json(joinpath(outdir,"summary.json"),summary)
    println("FINAL status=$status accepted=$(audit_result.verdict.accepted) reasons=$(audit_result.verdict.reasons)");flush(stdout)
    return summary
end
