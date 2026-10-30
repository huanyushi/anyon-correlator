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
        Dict("accepted"=>result.accepted,"status"=>result.status,"steps"=>result.steps)))
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
    evaluate=physical_evaluator(b,hx,hz;probe,chi)
    md["chi"]=chi;md["probe_sha256"]=bytes2hex(sha256(read(options["probe"])))
    summary=optimize_steps(x,b;evaluate,outdir=out,maxsteps=integer("steps",200),deadline,metadata=md)
    # A stationary candidate needs a fresh derivative-length and physical audit.
    final,_,_=load_state(joinpath(out,"final.json"))
    if summary.status=="stationary_candidate"
        return finalize_run(summary,final,b,md,out;deadline,assess=()->begin
            H,_=Legacy.operators(hx,hz);y=split_parameters(final,b)
            r=Legacy.dominant_contract(y,b.D,H;chi,tol=1e-10,maxiter=500,seeds=(1,17))
            _,g20=energy_gradient(final,b,H,r;method=:unrolled,unroll_steps=20)
            _,g40=energy_gradient(final,b,H,r;method=:unrolled,unroll_steps=40)
            stable=norm(g40-g20)<=2e-6 && norm(g40)<=1e-5
            stable || return (;status="gradient_unstable",accepted=false,reasons=["gradient_length_or_norm"],gradient_delta=norm(g40-g20))
            time()<deadline || return (;status="budget",accepted=false,reasons=["deadline"])
            audit=Legacy.dominant_audit(y,b.D,hx,hz;chi,chi_check=2chi,
                gradnorm=norm(g40),status="stationary",tol=1e-10,maxiter=500)
            return (;status=audit.verdict.accepted ? "stationary" : "audit_rejected",
                accepted=audit.verdict.accepted,reasons=audit.verdict.reasons,
                gradient_delta=norm(g40-g20),audit)
        end)
    end
    return summary
end
