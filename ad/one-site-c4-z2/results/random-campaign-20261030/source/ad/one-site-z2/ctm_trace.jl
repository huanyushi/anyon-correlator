# Diagnostic only: persist unconverged environments and their physical drift.
include("core.jl")
using .OneSiteZ2AD, SHA, Random
const A=OneSiteZ2AD
const L=A.Legacy
function main(args=ARGS)
    opts=Dict{String,String}(Tuple(split(replace(arg,r"^--"=>""),'=';limit=2)) for arg in args)
    path=opts["state"];out=opts["out"];chi=parse(Int,get(opts,"chi","16"));seed=parse(Int,get(opts,"seed","2026"))
    ispath(out) && error("diagnostic output exists")
    x,b,md=A.load_state(path);hx=Float64(md["hx"]);hz=Float64(md["hz"])
    mkpath(out)
    open(joinpath(out,"settings.json"),"w") do io
        A.JSON.print(io,(;input_sha256=bytes2hex(sha256(read(path))),chi,seed,hx,hz,D=b.D,tol=1e-10,maxiter=4000,block=100,accepted=false))
    end
    psi=L.make_peps(A.split_parameters(x,b),b.D);_,ops=L.operators(hx,hz)
    Random.seed!(seed)
    env=L.CTMRGEnv(randn,ComplexF64,psi,L.ComplexSpace(chi))
    for block in 1:40
        start=time()
        env,info=L.leading_boundary(env,psi,L.SimultaneousCTMRG(;tol=1e-10,maxiter=100,verbosity=0))
        m=L.measure(psi,env,ops,hx,hz)
        row=(;iteration_upper_bound=100block,converged=info.converged,residual=info.convergence_error,measurement=m,seconds=time()-start)
        open(joinpath(out,"trace.jsonl"),"a") do io;println(io,A.JSON.json(row));flush(io);end
        L.jldsave(joinpath(out,"env-$(100block).jld2");env,parameters=x,D=b.D,chi,seed,source_sha256=bytes2hex(sha256(read(path))))
        println("CTM TRACE upper_iter=$(100block) residual=$(info.convergence_error) e=$(m.energy) mx=$(m.mx) mz=$(m.mz) converged=$(info.converged)");flush(stdout)
        if info.converged
            open(joinpath(out,"summary.json"),"w") do io;A.JSON.print(io,row);end
            break
        end
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    main()
end
