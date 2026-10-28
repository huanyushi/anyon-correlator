using Zygote, Printf

struct ContractionFailure <: Exception
    message::String
end
Base.showerror(io::IO,e::ContractionFailure) = print(io,e.message)

function contract(x,D,H;chi=8,tol=1e-9,maxiter=300,seed=2026)
    psi = make_peps(x,D)
    Random.seed!(seed)
    env0 = CTMRGEnv(randn,ComplexF64,psi,ComplexSpace(chi))
    alg = SimultaneousCTMRG(;tol,maxiter,verbosity=0)
    start = time()
    env,info = leading_boundary(env0,psi,alg)
    info.converged && isfinite(info.convergence_error) || throw(ContractionFailure("CTMRG chi=$chi seed=$seed residual=$(info.convergence_error)"))
    value = expectation_value(psi,H,env)/8
    isfinite(real(value)) && abs(imag(value))<=1e-7 || throw(ContractionFailure("nonreal/nonfinite energy $value"))
    return (;energy=real(value),psi,env,info,seconds=time()-start,chi,seed)
end

function energy_gradient(x,D,H,result;tol=1e-9,maxiter=300,grad_tol=1e-8)
    alg = SimultaneousCTMRG(;tol,maxiter,verbosity=0)
    gad = PEPSKit.GradientAlgorithm(;alg=:FixedPointGradient,tol=grad_tol,maxiter=30,verbosity=1)
    # Two AD stages implement exactly the chain rule but avoid asking Zygote
    # to compile the entire CTMRG pullback and split-network graph together.
    value,grads = Zygote.withgradient(result.psi) do psi
        env,_ = PEPSKit.hook_pullback(leading_boundary,result.env,psi,alg;alg_rrule=gad)
        real(PEPSKit.cost_function(psi,env,H))/8
    end
    peps_gradient = only(grads)
    gradient = only(Zygote.gradient(x) do y
        psi = make_peps(y,D)
        real(sum(dot(peps_gradient.A[i],psi.A[i]) for i in eachindex(psi.A)))
    end)
    all(isfinite,gradient) || error("AD parameter gradient is nonfinite")
    abs(value-result.energy)<=1e-7 || throw(ContractionFailure("AD forward objective changed by $(value-result.energy)"))
    return Float64(value),Vector{Float64}(gradient)
end

function measure(psi,env,ops,hx,hz)
    values = Dict(name=>expectation_value(psi,op,env) for (name,op) in ops)
    imagmax = maximum(abs(imag(v)) for v in Base.values(values))
    vals = Dict(k=>Float64(real(v)) for (k,v) in values)
    astar = sum(vals["A$r$c"] for r in 1:2,c in 1:2)/4
    bplaq = sum(vals["B$r$c"] for r in 1:2,c in 1:2)/4
    mx = sum(vals["$k$r$c"] for k in ("XE","XN"),r in 1:2,c in 1:2)/8
    mz = sum(vals["$k$r$c"] for k in ("ZE","ZN"),r in 1:2,c in 1:2)/8
    return (;energy=-(astar+bplaq)/2-hx*mx-hz*mz,mx,mz,astar,bplaq,imagmax,observables=vals)
end

function physical_ok(m,hx,hz;tol=1e-6)
    all(isfinite,Base.values(m.observables)) && isfinite(m.energy) && m.imagmax<=1e-7 &&
        maximum(abs,Base.values(m.observables))<=1+tol && m.energy>=-1-abs(hx)-abs(hz)-tol
end

function calibration_ok(rows)
    length(rows)>=4 && all(row->row.status=="converged",rows) || return false
    return all(rows) do row
        m = row.measurement
        abs(m.energy+1)<=1e-8 && m.imagmax<=1e-7 &&
            maximum(abs(m.observables["$k$r$c"]-1) for k in ("A","B"),r in 1:2,c in 1:2)<=1e-8
    end
end

function probe_matches(p,D,chi,hx,hz,nparameters)
    get(p,"passed",false)==true && get(p,"D",nothing)==D &&
        get(p,"chi",nothing)==chi && get(p,"hx",nothing)==hx &&
        get(p,"hz",nothing)==hz && get(p,"parameter_count",nothing)==nparameters
end

function assess_audit(rows,hx,hz;gradnorm,status)
    reasons = String[]
    status=="stationary" && gradnorm<=1e-5 || push!(reasons,"optimizer_not_stationary")
    good = filter(r->r.status=="converged",rows)
    length(good)==length(rows) || push!(reasons,"contraction_failed")
    isempty(good) && return (;accepted=false,reasons=vcat(reasons,"no_contractions"))
    chi = minimum(r.chi for r in rows)
    base = filter(r->r.chi==chi,good)
    length(unique(r.seed for r in base))>=3 || push!(reasons,"insufficient_initializations")
    any(r->r.chi>chi,good) || push!(reasons,"missing_chi_check")
    all(r->physical_ok(r.measurement,hx,hz),good) || push!(reasons,"physical_bounds")
    measurements = [r.measurement for r in good]
    maximum(m.energy for m in measurements)-minimum(m.energy for m in measurements)<=1e-6 || push!(reasons,"energy_disagreement")
    keys0 = keys(first(measurements).observables)
    spread = maximum(maximum(m.observables[k] for m in measurements)-minimum(m.observables[k] for m in measurements) for k in keys0)
    spread<=1e-5 || push!(reasons,"observable_disagreement")
    ref = reference_values(hx,hz)
    all(m->abs(m.energy-ref.energy)<=1e-5,measurements) || push!(reasons,"energy_benchmark")
    for name in (:mx,:mz)
        target = getproperty(ref,name)
        limit = abs(target)>1e-12 ? 0.01abs(target) : 1e-6
        all(m->abs(getproperty(m,name)-target)<=limit,measurements) || push!(reasons,"$(name)_benchmark")
    end
    if hx==0 || hz==0
        prefix = hx==0 ? "B" : "A"
        all(m->maximum(abs(m.observables["$prefix$r$c"]-1) for r in 1:2,c in 1:2)<=1e-6,measurements) || push!(reasons,"conserved_sector")
    end
    if hx==0 && hz==0
        all(m->maximum(abs(m.observables["$k$r$c"]-1) for k in ("A","B"),r in 1:2,c in 1:2)<=1e-6,measurements) || push!(reasons,"zero_field_stabilizers")
    end
    return (;accepted=isempty(reasons),reasons)
end

function audit(x,D,hx,hz;chi=8,chi_check=16,tol=1e-9,maxiter=300,gradnorm=Inf,status="unknown")
    chi_check>chi || throw(ArgumentError("chi_check must exceed chi"))
    H,ops = operators(hx,hz)
    rows = NamedTuple[]
    for (ch,seed) in ((chi,2026),(chi,1),(chi,2),(chi_check,2026))
        row = try
            result = contract(x,D,H;chi=ch,tol,maxiter,seed)
            measurement = measure(result.psi,result.env,ops,hx,hz)
            (;chi=ch,seed,status="converged",residual=Float64(result.info.convergence_error),measurement)
        catch err
            err isa ContractionFailure || rethrow()
            (;chi=ch,seed,status="failed",error=sprint(showerror,err))
        end
        push!(rows,row)
        println("audit chi=$ch seed=$seed status=$(row.status)");flush(stdout)
    end
    return (;rows,verdict=assess_audit(rows,hx,hz;gradnorm,status))
end

function derivative_probe(x,D,hx,hz;chi=8,tol=1e-10,maxiter=300)
    H,ops = operators(hx,hz)
    t0 = time()
    result = contract(x,D,H;chi,tol,maxiter)
    compile_objective = time()-t0
    println("probe objective ready energy=$(result.energy) residual=$(result.info.convergence_error) wall=$(round(compile_objective;digits=2))s; starting two-stage AD");flush(stdout)
    t0 = time()
    energy,g = energy_gradient(x,D,H,result;tol,maxiter)
    compile_gradient = time()-t0
    println("probe first energy=$energy gradnorm=$(norm(g)) objective=$(round(compile_objective;digits=2))s gradient=$(round(compile_gradient;digits=2))s");flush(stdout)
    t0 = time(); repeat = contract(x,D,H;chi,tol,maxiter); objective_seconds = time()-t0
    t0 = time(); _,g2 = energy_gradient(x,D,H,repeat;tol,maxiter); gradient_seconds = time()-t0
    dirs = [normalize(randn(MersenneTwister(s),length(x))) for s in (17,31)]
    norm(g)>1e-12 && push!(dirs,-g/norm(g))
    checks = NamedTuple[]
    for (i,d) in enumerate(dirs),epsilon in (1e-4,1e-5)
        plus = contract(x+epsilon*d,D,H;chi,tol,maxiter).energy
        minus = contract(x-epsilon*d,D,H;chi,tol,maxiter).energy
        fd = (plus-minus)/(2epsilon)
        ad = dot(g,d)
        passed = abs(fd-ad)<=max(2e-5,0.005abs(ad))
        push!(checks,(;direction=i,epsilon,fd,ad,passed))
        println("FD direction=$i eps=$epsilon ad=$ad fd=$fd passed=$passed");flush(stdout)
    end
    return (;D,chi,hx,hz,parameter_count=length(x),energy,gradnorm=norm(g),measurement=measure(result.psi,result.env,ops,hx,hz),compile_objective,compile_gradient,objective_seconds,gradient_seconds,repeat_gradient_difference=norm(g-g2),checks,passed=all(c->c.passed,checks))
end
