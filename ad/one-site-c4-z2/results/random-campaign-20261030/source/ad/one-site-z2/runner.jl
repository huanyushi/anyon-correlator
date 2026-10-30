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
function run(options::AbstractDict)
    number(key,default)=parse(Float64,get(options,key,string(default)))
    integer(key,default)=parse(Int,get(options,key,string(default)))
    hx=number("hx",0.);hz=number("hz",0.1);chi=integer("chi",8)
    budget=number("budget",1800.);isfinite(budget) && budget>=0 || throw(ArgumentError("invalid budget"))
    deadline=min(time()+budget,parse(Float64,get(ENV,"AD_DEADLINE","Inf")))
    out=options["out"];probe=JSON.parsefile(options["probe"])
    if haskey(options,"resume")
        path=options["resume"];x,b,old=load_state(path);check_resume(old,hx,hz)
        haskey(options,"D") && integer("D",b.D)!=b.D && throw(ArgumentError("resume D mismatch"))
        md=merge(old,Dict("parent_sha256"=>bytes2hex(sha256(read(path))),"parent_path"=>abspath(path)))
    else
        b=ParameterBasis(integer("D",4));seed=integer("seed",71)
        init=get(options,"init","random");init in ("exact","random") || throw(ArgumentError("unknown initialization"))
        noise=number("noise",0.);isfinite(noise) && noise>=0 || throw(ArgumentError("invalid noise"))
        init=="random" && noise!=0 && throw(ArgumentError("noise is only for analytic initialization"))
        x=initial_parameters(b,seed;exact=init=="exact")
        noise>0 && (x+=noise*initial_parameters(b,seed))
        origin=init=="random" ? "symmetric_random" : noise>0 ? "perturbed_exact_zero_initialization" : "exact_zero_initialization"
        md=Dict{String,Any}("origin"=>origin,"seed"=>seed,"noise"=>noise,"hx"=>hx,"hz"=>hz)
    end
    tol=number("tol",1e-8);maxiter=integer("ctm-maxiter",300)
    energy_tolerance=number("energy-tolerance",1e-4)
    isfinite(tol) && 0<tol<=1e-6 || throw(ArgumentError("invalid CTMRG tolerance"))
    maxiter>0 && isfinite(energy_tolerance) && 0<energy_tolerance<=1e-3 || throw(ArgumentError("invalid accuracy settings"))
    evaluate=physical_evaluator(b,hx,hz;probe,chi,tol,maxiter)
    md["chi"]=chi;md["probe_sha256"]=bytes2hex(sha256(read(options["probe"])))
    md["ctm_tol"]=tol;md["ctm_maxiter"]=maxiter;md["energy_tolerance"]=energy_tolerance
    reference=Legacy.reference_values(hx,hz)
    summary=optimize_steps(x,b;evaluate,outdir=out,maxsteps=integer("steps",200),deadline,metadata=md,
        target_energy=reference.energy,energy_tolerance)
    # The energy benchmark is followed by a short independent contraction check.
    final,_,_=load_state(joinpath(out,"final.json"))
    if summary.status in ("stationary_candidate","energy_candidate")
        return finalize_run(summary,final,b,md,out;deadline,assess=()->energy_check(final,b,hx,hz;
            chi,chi_check=integer("chi-check",12),tol,maxiter,energy_tolerance,deadline,out))
    end
    return summary
end

"Energy-based finite-accuracy check; no stationarity or spectral claims."
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
