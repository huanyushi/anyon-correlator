"L-BFGS with an energy target, line search and per-update state saving."
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

# Run one initial PEPS at one field point; random.jl only schedules multiple seeds.
function run(options::AbstractDict)
    haskey(options,"resume") && throw(ArgumentError("use the analytic or random entry; --resume is not a public option"))
    number(key,default)=parse(Float64,get(options,key,string(default)))
    integer(key,default)=parse(Int,get(options,key,string(default)))
    hx=number("hx",0.);hz=number("hz",0.1);chi=integer("chi",8)
    budget=number("budget",1800.);isfinite(budget) && budget>=0 || throw(ArgumentError("invalid budget"))
    deadline=min(time()+budget,parse(Float64,get(ENV,"AD_DEADLINE","Inf")))
    out=options["out"]
    probe=haskey(options,"probe") ? JSON.parsefile(options["probe"]) : nothing
    b=ParameterBasis(integer("D",4));seed=integer("seed",71)
    init=get(options,"init","exact");init in ("exact","random") || throw(ArgumentError("unknown initialization"))
    noise=number("noise",init=="exact" ? 0.03 : 0.)
    isfinite(noise) && noise>=0 || throw(ArgumentError("invalid noise"))
    init=="random" && noise!=0 && throw(ArgumentError("noise is only for analytic initialization"))
    x=initial_parameters(b,seed;exact=init=="exact")
    noise>0 && (x+=noise*initial_parameters(b,seed))
    origin=init=="random" ? "symmetric_random" : noise>0 ? "perturbed_exact_zero_initialization" : "exact_zero_initialization"
    md=Dict{String,Any}("origin"=>origin,"seed"=>seed,"noise"=>noise,"hx"=>hx,"hz"=>hz)
    tol=number("tol",1e-8);maxiter=integer("ctm-maxiter",300)
    energy_tolerance=number("energy-tolerance",1e-4)
    isfinite(tol) && 0<tol<=1e-6 || throw(ArgumentError("invalid CTMRG tolerance"))
    maxiter>0 && isfinite(energy_tolerance) && 0<energy_tolerance<=1e-3 || throw(ArgumentError("invalid accuracy settings"))
    evaluate=physical_evaluator(b,hx,hz;probe,chi,tol,maxiter)
    md["chi"]=chi
    if probe!==nothing
        md["probe_sha256"]=bytes2hex(sha256(read(options["probe"])))
    end
    md["ctm_tol"]=tol;md["ctm_maxiter"]=maxiter;md["energy_tolerance"]=energy_tolerance
    reference=Legacy.reference_values(hx,hz)
    summary=optimize_steps(x,b;evaluate,outdir=out,maxsteps=integer("steps",200),deadline,metadata=md,
        target_energy=reference.energy,energy_tolerance)
    final,_,_=load_state(joinpath(out,"final.json"))
    if summary.status in ("stationary_candidate","energy_candidate")
        return finalize_run(summary,final,b,md,out;deadline,assess=()->energy_check(final,b,hx,hz;
            chi,chi_check=integer("chi-check",12),tol,maxiter,energy_tolerance,deadline,out))
    end
    return summary
end

function finalize_run(summary,x,b,md,out;deadline,assess)
    outcome=if time()>=deadline
        (;status="budget",accepted=false,reasons=["deadline"])
    else
        try
            assess()
        catch err
            err isa Legacy.ContractionFailure || err isa Legacy.GradientFailure || rethrow()
            (;status="audit_failed",accepted=false,reasons=[sprint(showerror,err)])
        end
    end
    time()>=deadline && (outcome=(;status="budget",accepted=false,reasons=["deadline"]))
    result=merge(summary,outcome)
    temporary=joinpath(out,"final-replacement.json")
    save_state(temporary,x,b;metadata=merge(Dict{String,Any}(string(k)=>v for (k,v) in pairs(md)),
        Dict("accepted"=>result.accepted,"status"=>result.status,"steps"=>result.steps,
             "energy"=>get(summary,:energy,nothing),"acceptance_kind"=>"finite_accuracy_series_energy")))
    mv(temporary,joinpath(out,"final.json");force=true)
    open(joinpath(out,"summary-replacement.json"),"w") do io
        JSON.print(io,result);flush(io)
    end
    mv(joinpath(out,"summary-replacement.json"),joinpath(out,"summary.json");force=true)
    return result
end

"Energy comparison and fresh CTMRG environments for the frozen candidate."
function energy_check(x,b,hx,hz;chi=8,chi_check=12,tol=1e-8,maxiter=300,
                      energy_tolerance=1e-4,deadline=Inf,out)
    chi_check>chi || throw(ArgumentError("check chi must exceed optimization chi"))
    H,ops=Legacy.operators(hx,hz);y=split_parameters(x,b)
    rows=NamedTuple[];reference=Legacy.reference_values(hx,hz)
    for ch in (chi,chi_check)
        time()<deadline || return (;status="budget",accepted=false,reasons=["deadline"])
        r=Legacy.dominant_contract(y,b.D,H;chi=ch,tol,maxiter,seed=2026,seeds=(71,),minimum_members=2)
        m=Legacy.measure(r.psi,r.env,ops,hx,hz)
        push!(rows,(;chi=ch,residual=r.info.convergence_error,measurement=m,dominance=r.dominance))
        open(joinpath(out,"energy-check-chi$ch.json"),"w") do io
            JSON.print(io,last(rows));flush(io)
        end
        println("ENERGY CHECK chi=$ch energy=$(m.energy) residual=$(r.info.convergence_error)");flush(stdout)
    end
    reasons=String[]
    all(r->Legacy.physical_ok(r.measurement,hx,hz),rows) || push!(reasons,"physical_bounds")
    all(r->abs(r.measurement.energy-reference.energy)<=energy_tolerance,rows) || push!(reasons,"series_energy")
    maximum(r.measurement.energy for r in rows)-minimum(r.measurement.energy for r in rows)<=1e-5 || push!(reasons,"chi_energy_difference")
    if hx==0
        all(r->maximum(abs(r.measurement.observables["B$i$j"]-1) for i in 1:2,j in 1:2)<=1e-3,rows) || push!(reasons,"conserved_B_sector")
    end
    accepted=isempty(reasons)
    return (;status=accepted ? "energy_accepted" : "energy_check_rejected",accepted,reasons,reference,
        energy_tolerance,acceptance_kind="finite_accuracy_series_energy",rows)
end
