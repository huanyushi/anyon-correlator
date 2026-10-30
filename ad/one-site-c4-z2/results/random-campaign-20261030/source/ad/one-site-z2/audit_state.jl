# Audit a frozen candidate, persisting every contraction and derivative stage.
module FrozenADAudit
include("core.jl")
using .OneSiteZ2AD, LinearAlgebra, SHA
const A=OneSiteZ2AD
function write_json(path,data)
    open(path*".tmp","w") do io;A.JSON.print(io,data);flush(io);end
    mv(path*".tmp",path;force=true)
end
function main(args=ARGS)
    opts=Dict{String,String}(Tuple(split(replace(arg,r"^--"=>""),'=';limit=2)) for arg in args)
    path=opts["state"];out=opts["out"]
    x,b,md=A.load_state(path);hx=Float64(md["hx"]);hz=Float64(md["hz"])
    digest=bytes2hex(sha256(read(path)))
    settings=(;input_sha256=digest,D=b.D,hx,hz,chis=[8,16],seeds=[2026,1,17],tol=1e-10,maxiter=1000)
    if ispath(out)
        get(opts,"resume-audit","false")=="true" || error("audit output exists")
        old=A.JSON.parsefile(joinpath(out,"settings.json"))
        old==A.JSON.parse(A.JSON.json(settings)) || error("audit settings or source changed")
    else
        mkpath(out);write_json(joinpath(out,"settings.json"),settings)
        A.save_state(joinpath(out,"input.json"),x,b;metadata=merge(md,Dict("accepted"=>false,"parent_sha256"=>digest)))
    end
    H,ops=A.Legacy.operators(hx,hz);y=A.split_parameters(x,b)
    gradient_path=joinpath(out,"gradient.json")
    for chi in (8,16),seed in (2026,1,17)
        rowpath=joinpath(out,"chi$(chi)-seed$(seed).json")
        if isfile(rowpath) && (chi!=8 || seed!=2026 || isfile(gradient_path));continue;end
        println("AUDIT START chi=$chi seed=$seed");flush(stdout)
        try
            r=A.contract(x,b,H;chi,seed,tol=1e-10,maxiter=1000)
            m=A.Legacy.measure(r.psi,r.env,ops,hx,hz)
            normcell=real(A.Legacy.network_value(r.psi,r.env))
            row=(;chi,seed,status="converged",residual=r.info.convergence_error,measurement=m,normcell)
            isfile(rowpath) || write_json(rowpath,row)
            println("AUDIT RESULT chi=$chi seed=$seed e=$(m.energy) mx=$(m.mx) mz=$(m.mz) A=$(m.astar) B=$(m.bplaq) residual=$(r.info.convergence_error)");flush(stdout)
            if chi==8 && seed==2026 && !isfile(gradient_path)
                gradients=[]
                for n in (20,40)
                    gp=joinpath(out,"gradient$n.json")
                    if isfile(gp)
                        push!(gradients,Float64.(A.JSON.parsefile(gp)["gradient"]))
                    else
                        value,g=A.energy_gradient(x,b,H,r;method=:unrolled,unroll_steps=n,tol=1e-10,maxiter=1000)
                        write_json(gp,(;value,gradient=g,norm=norm(g),unroll_steps=n));push!(gradients,g)
                    end
                    println("AUDIT GRADIENT n=$n norm=$(norm(last(gradients)))");flush(stdout)
                end
                delta=norm(gradients[1]-gradients[2]);gn=norm(gradients[2])
                write_json(gradient_path,(;delta,gradnorm=gn,stable=delta<=2e-6,stationary=gn<=1e-5))
            end
        catch err
            err isa A.Legacy.ContractionFailure || err isa A.Legacy.GradientFailure || rethrow()
            write_json(joinpath(out,"failure-chi$(chi)-seed$(seed).json"),(;error=sprint(showerror,err)))
            isfile(rowpath) || write_json(rowpath,(;chi,seed,status="failed",error=sprint(showerror,err)))
            println("AUDIT FAILED chi=$chi seed=$seed $(sprint(showerror,err))");flush(stdout)
        end
    end
    rows=map([(c,s) for c in (8,16) for s in (2026,1,17)]) do (chi,seed)
        d=A.JSON.parsefile(joinpath(out,"chi$(chi)-seed$(seed).json"))
        d["status"]=="converged" || return (;chi,seed,status="failed")
        m=d["measurement"]
        measurement=(;energy=m["energy"],mx=m["mx"],mz=m["mz"],astar=m["astar"],bplaq=m["bplaq"],imagmax=m["imagmax"],observables=m["observables"])
        (;chi,seed,status="converged",residual=d["residual"],normcell=d["normcell"],measurement)
    end
    g=isfile(gradient_path) ? A.JSON.parsefile(gradient_path) : Dict("stable"=>false,"stationary"=>false,"gradnorm"=>nothing)
    status=g["stable"] && g["stationary"] ? "stationary" : "gradient_failed"
    verdict=A.Legacy.assess_audit(rows,hx,hz;gradnorm=g["gradnorm"]===nothing ? Inf : g["gradnorm"],status)
    # Retain all independently initialized environments; never select by energy.
    norm_ok=all(c->begin
        ns=[r.normcell for r in rows if r.chi==c && r.status=="converged"]
        length(ns)==3 && all(>(0),ns) && (maximum(ns)-minimum(ns))/maximum(ns)<=1e-8
    end,(8,16))
    accepted=verdict.accepted && norm_ok
    summary=(;accepted,status=accepted ? "stationary" : "audit_rejected",verdict,norm_agreement=norm_ok,gradient=g,rows,input_sha256=digest)
    write_json(joinpath(out,"summary.json"),summary)
    if accepted && !isfile(joinpath(out,"accepted.json"))
        A.save_state(joinpath(out,"accepted.json"),x,b;metadata=merge(md,Dict("accepted"=>true,"status"=>"stationary",
            "parent_sha256"=>digest,"audit_sha256"=>bytes2hex(sha256(read(joinpath(out,"summary.json")))))))
    end
    println("AUDIT FINAL accepted=$accepted reasons=$(verdict.reasons) norm_agreement=$norm_ok");flush(stdout)
end
end
if abspath(PROGRAM_FILE)==@__FILE__
    FrozenADAudit.main()
end
