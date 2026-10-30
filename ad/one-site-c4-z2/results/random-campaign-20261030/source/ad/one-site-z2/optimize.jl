"Budgeted L-BFGS mechanics; physical acceptance is a separate audit."
function optimize_steps(initial,b::ParameterBasis;evaluate,outdir,maxsteps=200,
                         deadline=time()+1800,metadata=(;),target_energy=nothing,energy_tolerance=1e-4)
    ispath(outdir) && throw(ArgumentError("output already exists"))
    maxsteps>=0 || throw(ArgumentError("negative step count"))
    x=copy(initial);tensors(x,b);mkpath(outdir)
    md=Dict{String,Any}(string(k)=>v for (k,v) in pairs(metadata))
    md["accepted"]=false
    save_state(joinpath(outdir,"initial.json"),x,b;metadata=md)
    history=Tuple{Vector{Float64},Vector{Float64}}[]
    previous=nothing;steps=0;status="step_limit";start=time();last_energy=nothing
    record(row)=Legacy.write_step(outdir,row)
    try
        while steps<maxsteps
            if time()>=deadline;status="budget";break;end
            value=evaluate(x;gradient=true)
            last_energy=value.energy
            if time()>=deadline;status="budget";break;end
            g=value.gradient
            length(g)==length(x) && all(isfinite,g) && isfinite(value.energy) ||
                throw(Legacy.GradientFailure("invalid objective or gradient"))
            record((;step=steps,energy=value.energy,gradnorm=norm(g),seconds=time()-start))
            println("step=$steps energy=$(value.energy) gradnorm=$(norm(g))");flush(stdout)
            if target_energy!==nothing && abs(value.energy-target_energy)<=energy_tolerance
                status="energy_candidate";break
            end
            if norm(g)<=1e-5;status="stationary_candidate";break;end
            if previous!==nothing
                push!(history,(x-previous.x,g-previous.g))
                length(history)>8 && popfirst!(history)
            end
            d=Legacy.lbfgs_direction(g,history);alpha=min(1.,0.2/norm(d));accepted=false
            for attempt in 1:16
                if time()>=deadline;status="budget";break;end
                trial=x+alpha*d
                evaluation=try
                    checked_trial(trial,b)
                    evaluate(trial;gradient=false)
                catch err
                    err isa Legacy.ContractionFailure || rethrow()
                    record((;step=steps,attempt,event="contraction_rejected",message=sprint(showerror,err)))
                    nothing
                end
                if time()>=deadline;status="budget";break;end
                if evaluation!==nothing && isfinite(evaluation.energy) &&
                    evaluation.energy<=value.energy+1e-4alpha*dot(g,d)
                    previous=(;x=copy(x),g=copy(g));x=trial;steps+=1;last_energy=evaluation.energy
                    save_state(joinpath(outdir,"step$(lpad(steps,4,'0')).json"),x,b;
                        metadata=merge(md,Dict("step"=>steps,"energy"=>evaluation.energy)))
                    accepted=true;break
                end
                alpha/=2
            end
            if !accepted
                status=="budget" || (status="line_search_failed")
                break
            end
        end
    catch err
        err isa Legacy.GradientFailure || err isa Legacy.ContractionFailure || rethrow()
        status=err isa Legacy.GradientFailure ? "gradient_failed" : "contraction_failed"
        record((;step=steps,event=status,message=sprint(showerror,err)))
    end
    time()>=deadline && (status="budget")
    save_state(joinpath(outdir,"final.json"),x,b;metadata=merge(md,Dict("status"=>status,"steps"=>steps,"energy"=>last_energy)))
    summary=(;status,steps,energy=last_energy,accepted=false,seconds=time()-start)
    open(joinpath(outdir,"summary.json"),"w") do io
        JSON.print(io,summary);flush(io)
    end
    return summary
end
