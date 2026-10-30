# One persistent Julia process: compilation is shared across random starts.
include("core.jl")
using .OneSiteZ2AD, SHA, Dates
const A=OneSiteZ2AD

function campaign(args=ARGS)
    opts=Dict{String,String}(Tuple(split(arg[3:end],'=';limit=2)) for arg in args)
    out=opts["out"];ispath(out) && error("campaign output already exists")
    mkpath(out)
    budget=parse(Float64,get(opts,"budget","21600"))
    deadline=min(time()+budget,parse(Float64,get(ENV,"AD_DEADLINE","Inf")))
    seedlist=parse.(Int,split(get(opts,"seeds","23,41,67,89,113,137"),','))
    isempty(seedlist) && error("at least one fresh seed is required")
    resume=get(opts,"resume",joinpath(@__DIR__,"results/d4-random19-cont1-20261029/step0003.json"))
    probe=get(opts,"probe",joinpath(@__DIR__,"results/d4-gradient-probe-20261029/summary.json"))
    settings=(;hx=0.,hz=0.1,D=4,chi=8,chi_check=12,tol=1e-8,maxiter=300,
        energy_tolerance=1e-4,budget,seeds=seedlist,resume=abspath(resume),probe=abspath(probe),
        parent_sha256=bytes2hex(sha256(read(resume))),probe_sha256=bytes2hex(sha256(read(probe))),
        started=string(now()),initial_segment=3600,later_segment=2700)
    writejson(path,data)=open(path*".tmp","w") do io
        A.JSON.print(io,data);flush(io)
        # Rename after closing below, through the wrapper.
    end
    function persist(name,data)
        path=joinpath(out,name);writejson(path,data);mv(path*".tmp",path;force=true)
    end
    persist("settings.json",settings)
    snapshot=joinpath(out,"source");mkpath(snapshot)
    hashes=Dict{String,String}()
    for dir in (@__DIR__,joinpath(@__DIR__,"../split"),joinpath(@__DIR__,"../../model"))
        for file in readdir(dir;join=true)
            isfile(file) && (endswith(file,".jl") || endswith(file,".py")) || continue
            relative=relpath(file,joinpath(@__DIR__,"../.."))
            dest=joinpath(snapshot,relative);mkpath(dirname(dest));cp(file,dest)
            hashes[relative]=bytes2hex(sha256(read(file)))
        end
    end
    persist("source-hashes.json",hashes)
    attempts=Any[];best=nothing;nextseed=1;continuations=Dict{Int,Int}()
    start=time();attempt=0;accepted=false
    while time()<deadline-120
        attempt+=1
        # First finish the existing random-origin trajectory. Thereafter alternate
        # independent starts with the best trajectory, but only when it progressed.
        continuation=attempt==1 || (iseven(attempt-1) && best!==nothing &&
            best.improvement>=1e-3 && get(continuations,best.seed,0)<2)
        if continuation
            path=attempt==1 ? resume : best.path
            _,_,md=A.load_state(path);seed=Int(md["seed"])
            continuations[seed]=get(continuations,seed,0)+1
        else
            nextseed>length(seedlist) && break
            seed=seedlist[nextseed];nextseed+=1;path=nothing
        end
        segment=min(attempt==1 ? 3600. : 2700.,deadline-time()-30)
        dir=joinpath(out,"attempt$(lpad(attempt,2,'0'))-seed$seed")
        runopts=Dict("D"=>"4","chi"=>"8","chi-check"=>"12","hx"=>"0","hz"=>"0.1",
            "tol"=>"1e-8","ctm-maxiter"=>"300","energy-tolerance"=>"1e-4",
            "steps"=>"1000","budget"=>string(segment),"out"=>dir,"probe"=>probe)
        if path===nothing
            runopts["init"]="random";runopts["seed"]=string(seed)
        else
            runopts["resume"]=path
        end
        persist("progress.json",(;status="running",attempt,seed,continuation,active_directory=dir,
            attempts,best,seconds=time()-start,deadline))
        println("CAMPAIGN attempt=$attempt seed=$seed continuation=$continuation budget=$segment");flush(stdout)
        result=try
            A.run(runopts)
        catch err
            # LAPACK and programming failures remain explicit; don't mistake them
            # for a converged trajectory or destroy other completed attempts.
            message=sprint(showerror,err,catch_backtrace())
            persist("failure-attempt$attempt.json",(;message))
            println("CAMPAIGN attempt failed: $message");flush(stdout)
            (;status="worker_error",accepted=false)
        end
        rows=isfile(joinpath(dir,"trajectory.jsonl")) ?
            [A.JSON.parse(line) for line in eachline(joinpath(dir,"trajectory.jsonl"))] : Any[]
        energies=[Float64(r["energy"]) for r in rows if haskey(r,"gradnorm") && haskey(r,"energy")]
        checkpoints=isdir(dir) ? sort(filter(p->occursin(r"step\d{4}\.json$",p),readdir(dir;join=true))) : String[]
        final=isfile(joinpath(dir,"final.json")) ? joinpath(dir,"final.json") :
            !isempty(checkpoints) ? last(checkpoints) : isfile(joinpath(dir,"initial.json")) ? joinpath(dir,"initial.json") : nothing
        # Checkpoint metadata can contain a newer accepted update than the last gradient row.
        energy=isempty(energies) ? nothing : last(energies)
        if final!==nothing
            _,_,md=A.load_state(final)
            checkpoint_energy=get(md,"energy",nothing)
            checkpoint_energy===nothing || (energy=Float64(checkpoint_energy))
        end
        improvement=isempty(energies) || energy===nothing ? 0. : first(energies)-energy
        row=(;attempt,seed,continuation,status=result.status,accepted=result.accepted,
            energy,improvement,updates=length(checkpoints),path=final)
        push!(attempts,row)
        if energy!==nothing && final!==nothing && (best===nothing || energy<best.energy)
            best=row
        end
        accepted=result.accepted
        persist("progress.json",(;status=accepted ? "energy_accepted" : "between_attempts",
            attempts,best,seconds=time()-start,deadline))
        println("CAMPAIGN RESULT $(A.JSON.json(row))");flush(stdout)
        accepted && break
        GC.gc()
    end
    persist("summary.json",(;status=accepted ? "energy_accepted" : "budget_or_attempt_limit",
        accepted,attempts,best,seconds=time()-start,reference=A.Legacy.reference_values(0.,0.1)))
end

if abspath(PROGRAM_FILE)==@__FILE__
    campaign()
end
