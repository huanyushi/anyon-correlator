module TransferRun
const PROCESS_START=time()
include("core.jl")
using .ToricTransferSpectrum, SHA
const T=ToricTransferSpectrum
const S=T.OneSiteZ2SU

function options(args)
    d=Dict("state"=>"","out"=>"","candidate"=>"false","chi"=>"8","compare-chi"=>"12",
        "nk"=>"64","num"=>"8","budget"=>"7080","tol"=>"1e-6","fixedpoint-tol"=>"1e-5",
        "spectrum-tol"=>"1e-10","residual-tol"=>"1e-7","resolution"=>"1e-12")
    for arg in args
        kv=split(replace(arg,r"^--"=>""),'=';limit=2)
        length(kv)==2 && haskey(d,kv[1]) || throw(ArgumentError("unknown option $arg"))
        d[kv[1]]=kv[2]
    end
    all(k->!isempty(d[k]),("state","out")) || throw(ArgumentError("--state and --out are required"))
    d["candidate"] in ("true","false") || throw(ArgumentError("candidate must be true or false"))
    for key in ("budget","tol","fixedpoint-tol","spectrum-tol","residual-tol","resolution")
        v=parse(Float64,d[key]);isfinite(v) && v>0 || throw(ArgumentError("invalid $key"))
    end
    for key in ("nk","num")
        parse(Int,d[key])>0 || throw(ArgumentError("$key must be positive"))
    end
    chi=parse(Int,d["chi"]);compare=parse(Int,d["compare-chi"])
    T.boundary_chis((;chis=(chi,)))
    compare==0 || (T.boundary_chis((;chis=(compare,)));compare>chi) ||
        throw(ArgumentError("compare-chi must be zero or larger than chi"))
    return d
end

function main(args=ARGS)
    opts=options(args);out=opts["out"]
    deadline=PROCESS_START+parse(Float64,opts["budget"])
    chi=parse(Int,opts["chi"]);checkchi=parse(Int,opts["compare-chi"])
    nk=parse(Int,opts["nk"]);num=parse(Int,opts["num"])
    tol=parse(Float64,opts["tol"]);fixedpoint_tol=parse(Float64,opts["fixedpoint-tol"])
    spectrum_tol=parse(Float64,opts["spectrum-tol"]);residual_tol=parse(Float64,opts["residual-tol"])
    resolution=parse(Float64,opts["resolution"])
    candidate=opts["candidate"]=="true"
    st,metadata=T.load_production_state(opts["state"];reference=candidate)
    S.prepare_output(out)
    writejson(name,data)=S.write_json(joinpath(out,name),data)
    cp(opts["state"],joinpath(out,"input.json"))
    source=joinpath(out,"source");mkpath(source)
    hashes=Dict{String,String}()
    for name in ("core.jl","numerics.jl","input.jl","network.jl","boundary.jl","spectrum.jl","output.jl","run.jl","plot.py")
        path=joinpath(@__DIR__,name);cp(path,joinpath(source,name))
        hashes[name]=bytes2hex(sha256(read(path)))
    end
    ks=collect(range(-pi;step=2pi/nk,length=nk))
    writejson("settings.json",(;opts,input_sha256=bytes2hex(sha256(read(opts["state"]))),input_metadata=metadata,
        candidate_input=candidate,direction="down",step=1,momenta=ks,num,chi,checkchi,
        boundary_tol=tol,fixedpoint_tol,spectrum_tol,residual_tol,resolution,source_sha256=hashes,
        julia=string(VERSION),mpskit=string(pkgversion(T.MPSKit)),
        scope="frozen analytic-origin optimized PEPS; qualitative magnitude spectrum; not Hamiltonian excitations"))
    jsonl=joinpath(out,"spectrum.jsonl");touch(jsonl)
    points=SectorPoint[];boundary=nothing;grid_complete=false;comparison=nothing
    function summarize(status;message=nothing)
        residuals=[p.projected_residual for p in points if p.projected_residual!==nothing]
        writejson("summary.json",(;status,grid_complete,count=length(points),expected=nk*4num,
            candidate_input=candidate,input_ground_state_accepted=get(metadata,"accepted",false),
            max_projected_residual=isempty(residuals) ? nothing : maximum(residuals),
            chi,comparison,message,boundary,seconds=time()-PROCESS_START,peak_rss_bytes=Sys.maxrss(),
            interpretation="finite-chi single-site quasiparticle variational transfer spectrum"))
    end
    solve_options=(;tol=spectrum_tol,residual_tol,resolution,deadline)
    function point_batch(net,pair,k;path=jsonl,destination=points)
        batch=SectorPoint[]
        for g in 0:1,q in 0:1
            sector=Base.invokelatest(solve_sector,net,pair,SectorLabel(g,q);momenta=[k],num,options=solve_options)
            # Preserve completed sectors even if the next one reaches the deadline.
            for p in sector
                T.append_point(path,p);push!(destination,p);push!(batch,p)
            end
        end
        return batch
    end
    # Solver precision and finite-chi accuracy are recorded separately.
    point_ok(p)=p.status in ("converged","below_resolution","zero_mode") &&
        p.projected_residual!==nothing && isfinite(p.projected_residual) && p.projected_residual<=residual_tol
    try
        println("TRANSFER chi=$chi nk=$nk modes=$num candidate=$candidate");flush(stdout)
        net=Base.invokelatest(build_transfer,st)
        pair=Base.invokelatest(solve_boundaries,net;options=(;chis=(chi,),tol,fixedpoint_tol,maxiter=300,deadline))
        T.save_boundaries(joinpath(out,"boundaries.json"),pair);boundary=pair.diagnostics
        for (index,k) in enumerate(ks)
            time()<deadline || error("main spectrum time budget reached")
            started=time();batch=point_batch(net,pair,k)
            T.write_csv(jsonl,joinpath(out,"spectrum.csv"))
            writejson("progress.json",(;stage="grid",momentum=index,total=nk,k,count=length(points),seconds=time()-PROCESS_START))
            println("GRID $index/$nk k/pi=$(k/pi) points=$(length(batch)) seconds=$(time()-started)");flush(stdout)
            all(point_ok,batch) || error("a grid point failed its projected residual check")
        end
        grid_complete=length(points)==nk*4num
        summarize("grid_finished")
        # Optional chi comparison never blocks saving a complete main grid.
        if checkchi>0 && deadline-time()>300
            compare_rows=SectorPoint[];deltas=NamedTuple[]
            try
                writejson("progress.json",(;stage="chi_comparison",chi=checkchi,seconds=time()-PROCESS_START))
                checked=Base.invokelatest(solve_boundaries,net;options=(;chis=(checkchi,),tol,fixedpoint_tol,maxiter=300,deadline))
                T.save_boundaries(joinpath(out,"boundaries-chi$checkchi.json"),checked)
                comparefile=joinpath(out,"comparison-chi$checkchi.jsonl");touch(comparefile)
                for k in (0.,pi/2,-pi)
                    # Match the closest actual grid point, without relying on Float64 equality.
                    index=argmin(abs.(ks.-k));abs(ks[index]-k)<1e-12 || continue
                    batch=point_batch(net,checked,ks[index];path=comparefile,destination=compare_rows)
                    T.write_csv(comparefile,joinpath(out,"comparison-chi$checkchi.csv"))
                    all(point_ok,batch) || error("chi comparison projected residual check failed")
                    for g in 0:1,q in 0:1
                        a=sort(filter(p->p.g==g && p.q==q && abs(p.k-ks[index])<1e-12,points);by=p->p.branch)
                        b=sort(filter(p->p.g==g && p.q==q,batch);by=p->p.branch)
                        differences=[abs(log(x.rho/y.rho)) for (x,y) in zip(a,b) if min(x.rho,y.rho)>resolution]
                        push!(deltas,(;k=ks[index],g,q,leading=isempty(differences) ? nothing : first(differences),
                            maximum=isempty(differences) ? nothing : maximum(differences),compared=length(differences)))
                    end
                    comparison=(;status="partial",chi=checkchi,count=length(compare_rows),deltas)
                    writejson("chi-comparison.json",comparison);summarize("grid_finished")
                    println("CHI COMPARISON chi=$checkchi k/pi=$(k/pi) count=$(length(compare_rows))");flush(stdout)
                end
                comparison=(;status="finished",chi=checkchi,count=length(compare_rows),deltas)
            catch err
                if @isdefined(comparefile) && isfile(comparefile)
                    T.write_csv(comparefile,joinpath(out,"comparison-chi$checkchi.csv"))
                end
                comparison=(;status=time()>=deadline ? "budget" : "error",chi=checkchi,count=length(compare_rows),deltas,error=sprint(showerror,err))
            end
        else
            comparison=(;status="skipped",reason="disabled or fewer than 300 seconds remaining")
        end
        writejson("chi-comparison.json",comparison)
        summarize("finished")
    catch err
        T.write_csv(jsonl,joinpath(out,"spectrum.csv"))
        summarize(time()>=deadline ? "budget" : "error";message=sprint(showerror,err))
        rethrow()
    end
    println("TRANSFER END grid_complete=$grid_complete count=$(length(points))");flush(stdout)
end
end
if abspath(PROGRAM_FILE)==@__FILE__
    TransferRun.main()
end
