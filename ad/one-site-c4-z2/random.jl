# Seed scheduling only. Every initial state uses the shared optimizer in core.jl.
include("core.jl")

function random_search(options::AbstractDict)
    A=OneSiteZ2AD
    seeds=parse.(Int,split(get(options,"seeds",get(options,"seed","19")),','))
    isempty(seeds) && throw(ArgumentError("at least one random seed is required"))
    length(unique(seeds))==length(seeds) || throw(ArgumentError("duplicate random seeds"))
    # A single seed is exactly one ordinary optimization, with the whole budget.
    if length(seeds)==1
        return A.run(merge(options,Dict("init"=>"random","seed"=>string(only(seeds)))))
    end
    budget=parse(Float64,get(options,"budget","1800"))
    isfinite(budget) && budget>0 || throw(ArgumentError("invalid budget"))
    out=options["out"];ispath(out) && throw(ArgumentError("output already exists"))
    deadline=min(time()+budget,parse(Float64,get(ENV,"AD_DEADLINE","Inf")))
    mkpath(out)
    function persist(name,data)
        path=joinpath(out,name)
        open(path*".tmp","w") do io
            A.JSON.print(io,data);flush(io)
        end
        mv(path*".tmp",path;force=true)
    end
    hx=parse(Float64,get(options,"hx","0"));hz=parse(Float64,get(options,"hz","0.1"))
    persist("settings.json",(;seeds,budget,hx,hz,D=parse(Int,get(options,"D","4")),
        allocation="remaining time divided between remaining seeds"))
    attempts=Any[];best=nothing;accepted=false;start=time()
    for (index,seed) in enumerate(seeds)
        remaining=deadline-time()
        remaining>0 || break
        allocation=remaining/(length(seeds)-index+1)
        dir=joinpath(out,"attempt$(lpad(index,2,'0'))-seed$seed")
        current=merge(options,Dict("init"=>"random","seed"=>string(seed),
            "budget"=>string(allocation),"out"=>dir))
        persist("progress.json",(;status="running",seed,active_directory=dir,attempts,best))
        println("RANDOM seed=$seed budget=$allocation");flush(stdout)
        result=A.run(current)
        energy=result.energy
        row=(;seed,status=result.status,accepted=result.accepted,energy,
            steps=result.steps,path=joinpath(dir,"final.json"))
        push!(attempts,row)
        if energy!==nothing && (best===nothing || energy<best.energy)
            best=row
        end
        accepted=result.accepted
        persist("summary.json",(;status=accepted ? "energy_accepted" : "incomplete",
            accepted,attempts,best,seconds=time()-start))
        persist("progress.json",(;status="between_seeds",attempts,best))
        accepted && break
        GC.gc()
    end
    summary=(;status=accepted ? "energy_accepted" : "finished_without_accepted_state",
        accepted,attempts,best,seconds=time()-start)
    persist("summary.json",summary)
    return summary
end

if abspath(PROGRAM_FILE)==@__FILE__
    options=Dict{String,String}()
    for arg in ARGS
        startswith(arg,"--") && occursin('=',arg) || error("expected --key=value")
        key,value=split(arg[3:end],'=';limit=2);options[key]=value
    end
    get(options,"init","random")=="random" || error("analytic initial states use run.jl")
    random_search(options)
end
