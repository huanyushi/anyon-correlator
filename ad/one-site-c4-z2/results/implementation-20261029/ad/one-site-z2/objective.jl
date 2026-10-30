function checked_legacy_gradient(f)
    try
        f()
    catch err
        err isa ErrorException && err.msg=="AD parameter gradient is nonfinite" &&
            throw(Legacy.GradientFailure(err.msg))
        rethrow()
    end
end
function checked_trial(x,b)
    try
        tensors(x,b)
    catch err
        err isa ArgumentError && err.msg in ("unsupported weight range","zero tensor norm","invalid AD parameters") &&
            throw(Legacy.ContractionFailure(err.msg))
        rethrow()
    end
end
contract(x,b::ParameterBasis,H;kwargs...)=Legacy.contract(split_parameters(x,b),b.D,H;kwargs...)
function energy_gradient(x,b::ParameterBasis,H,reference;kwargs...)
    value,g=checked_legacy_gradient(()->Legacy.energy_gradient(split_parameters(x,b),b.D,H,reference;kwargs...))
    return value,parameter_pullback(x,b,g)
end
unrolled_energy(x,b::ParameterBasis,H,reference;kwargs...)=
    Legacy.unrolled_energy(split_parameters(x,b),b.D,H,reference;kwargs...)

function probe_matches(p,b,hx,hz;chi=8,unroll_steps=20)
    get(p,"charges",nothing)==b.charges &&
        Legacy.probe_matches(p,b.D,chi,hx,hz,parameter_count(b);
            method=:unrolled,unroll_steps)
end
function check_resume(md,hx,hz)
    get(md,"hx",nothing)==hx && get(md,"hz",nothing)==hz &&
        haskey(md,"origin") && haskey(md,"seed") ||
        throw(ArgumentError("resume requires matching fields and initialization provenance"))
    return md
end
function physical_evaluator(b,hx,hz;probe,chi=8,unroll_steps=20,tol=1e-10,maxiter=500)
    probe_matches(probe,b,hx,hz;chi,unroll_steps) || throw(ArgumentError("gradient evidence does not match setup"))
    H,ops=Legacy.operators(hx,hz)
    function evaluate(x;gradient=true)
        y=split_parameters(x,b)
        r=Legacy.dominant_contract(y,b.D,H;chi,tol,maxiter,seeds=(1,17),minimum_members=2)
        m=Legacy.measure(r.psi,r.env,ops,hx,hz)
        Legacy.physical_ok(m,hx,hz) && abs(m.energy-r.energy)<=1e-8 ||
            throw(Legacy.ContractionFailure("invalid observables or energy normalization"))
        if gradient
            value,g=energy_gradient(x,b,H,r;method=:unrolled,unroll_steps,tol,maxiter)
            return (;energy=value,gradient=g,measurement=m,dominance=r.dominance)
        end
        return (;energy=r.energy,gradient=nothing,measurement=m,dominance=r.dominance)
    end
    return evaluate
end
