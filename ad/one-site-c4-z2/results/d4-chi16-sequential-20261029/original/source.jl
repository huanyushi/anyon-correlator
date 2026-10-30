include("/home/huanyushi/code/anyon-correlator/ad/one-site-z2/core.jl")
using .OneSiteZ2AD, PEPSKit, JLD2, SHA
const A=OneSiteZ2AD
function main()
    root="/home/huanyushi/code/anyon-correlator/"
    source=root*"ad/one-site-z2/results/d4-exact-noise-cont1-20261029/final.json"
    checkpoint=root*"ad/one-site-z2/results/d4-chi16-trace-20261029/env-4000.jld2"
    out=root*"ad/one-site-z2/results/d4-chi16-sequential-20261029"
    ispath(out) && error("output exists")
    mkpath(out)
    cp(@__FILE__,joinpath(out,"source.jl"))
    x,b,md=A.load_state(source)
    psi=A.Legacy.make_peps(A.split_parameters(x,b),b.D)
    _,ops=A.Legacy.operators(md["hx"],md["hz"])
    env=JLD2.load(checkpoint,"env")
    open(joinpath(out,"settings.json"),"w") do io
        A.JSON.print(io,(;algorithm="SequentialCTMRG",tol=1e-10,chi=16,maxiter=300,
            input_sha256=bytes2hex(sha256(read(source))),environment_sha256=bytes2hex(sha256(read(checkpoint))),accepted=false))
    end
    for block in 1:6
        env,info=leading_boundary(env,psi,SequentialCTMRG(;tol=1e-10,maxiter=50,verbosity=0))
        measurement=A.Legacy.measure(psi,env,ops,md["hx"],md["hz"])
        row=(;iteration_upper_bound=50block,residual=info.convergence_error,converged=info.converged,measurement)
        open(joinpath(out,"trace.jsonl"),"a") do io;println(io,A.JSON.json(row));end
        jldsave(joinpath(out,"env-$(50block).jld2");env)
        println(A.JSON.json(row));flush(stdout)
        info.converged && break
    end
end
main()
