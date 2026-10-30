module FilteredFigure7
const START=time()
include("../../transfer-spectrum/numerics.jl")
using JSON
include("../filtered-toric-code/tensor.jl")
phase_flip(r)=InfiniteMPS([-r.AL[1]],[-r.AR[1]],[r.C[1]],[-r.AC[1]])
include("topological.jl")
include("charge_condensed.jl")
include("flux_condensed.jl")
const PANELS=Dict("a"=>(.55,.45,:topological),"b"=>(.8,.45,:topological),
    "c"=>(1.05,.45,:charge_condensed),"d"=>(1.3,.45,:charge_condensed),
    "e"=>(.45,.55,:topological),"f"=>(.45,.8,:topological),
    "g"=>(.45,1.05,:flux_condensed),"h"=>(.45,1.3,:flux_condensed))
phase_setup(phase,A,chi;kwargs...)=getfield(@__MODULE__,Symbol(phase,"_setup"))(A,chi;kwargs...)
write_json(path,data)=open(io->(JSON.print(io,data);println(io)),path,"w")

# Same full QP matrix as Fig.7(a), accepting Trivial, Z2 or Z2×Z2 charge.
function sector_modes(op,l,r,sector,k;num,tol,deadline)
    proto=Base.invokelatest(MPSKit.LeftGaugedQP,randn,l,r;sector,momentum=Float64(k))
    n=sum(dim,proto.Xs)
    action,_,_=Base.invokelatest(effective_action,op,l,r,proto;tol)
    matrix=Matrix{ComplexF64}(undef,n,n)
    for j in 1:n
        time()<deadline || error("spectrum time budget reached")
        x=zeros(ComplexF64,n);x[j]=1
        matrix[:,j]=qp_coordinates(action(qp_from_coordinates(proto,x)))
    end
    eig=eigen(matrix);order=sortperm(abs.(eig.values);rev=true)[1:min(num,n)]
    values=eig.values[order]
    residuals=[norm(action(qp_from_coordinates(proto,eig.vectors[:,j]))-
        eig.values[j]*qp_from_coordinates(proto,eig.vectors[:,j])) for j in order]
    (;values,residuals,dimension=n)
end

function main(args=ARGS)
    opts=Dict("panel"=>"b","out"=>"","chi"=>"16","nk"=>"41","num"=>"64",
        "tol"=>"1e-6","spectrum-tol"=>"1e-7","budget"=>"1800")
    for arg in args
        k,v=split(replace(arg,r"^--"=>""),'=';limit=2)
        haskey(opts,k) || error("unknown option $k");opts[k]=v
    end
    panel=opts["panel"];bz,bx,phase=PANELS[panel];out=opts["out"]
    isempty(out) && error("--out required");ispath(out) && error("output exists")
    chi=parse(Int,opts["chi"]);nk=parse(Int,opts["nk"]);num=parse(Int,opts["num"])
    tol=parse(Float64,opts["tol"]);stol=parse(Float64,opts["spectrum-tol"])
    deadline=START+parse(Float64,opts["budget"]);mkpath(out)
    root=normpath(joinpath(@__DIR__,"../.."))
    sources=vcat(["transfer-spectrum/$f.jl" for f in ("numerics","network","boundary","spectrum")],
        ["benchmarks/filtered-toric-code/tensor.jl"],
        ["benchmarks/filtered-toric-code-fig7/$f" for f in ("run.jl","topological.jl","charge_condensed.jl","flux_condensed.jl","plot.py")],
        ["julia-env/Project.toml","julia-env/Manifest.toml"])
    for file in sources
        dest=joinpath(out,"source",file);mkpath(dirname(dest));cp(joinpath(root,file),dest)
    end
    ks=collect(range(-pi,pi;length=nk))
    write_json(joinpath(out,"settings.json"),(;panel,betaz=bz,betax=bx,phase,chi,nk,num,tol,spectrum_tol=stol,momenta=ks,opts))
    count=0;maxres=0.
    try
        A=filtered_tensor(bx,bz)
        write_json(joinpath(out,"tensor.json"),(;shape=size(A),A=vec(A)))
        println("BOUNDARY panel=$panel phase=$phase chi=$chi");flush(stdout)
        s=phase_setup(phase,A,chi;tol,deadline)
        overlaps=[boundary_overlap(s.r,apply_group(s.r,g)) for g in (s.actions[1],s.actions[2],prod(s.actions))]
        a=convert(Array,s.r.AL[1]);d=density(s.r,s.op,s.env)
        write_json(joinpath(out,"boundary.json"),(;shape=size(a),AL_real=vec(real.(a)),AL_imag=vec(imag.(a)),
            virtual_space=string(space(s.r.AL[1],1)),basis=s.basis,residual=s.res,
            density_real=real(d),density_imag=imag(d),overlaps_ket_bra_both=overlaps))
        println("BOUNDARY residual=$(s.res) overlaps=$overlaps");flush(stdout)
        write_json(joinpath(out,"sectors.json"),[(;name=c.name,color=c.color,marker=c.marker,charge=string(c.sector)) for c in s.cases])
        open(joinpath(out,"spectrum.csv"),"w") do io
            println(io,"sector,k,branch,ratio_real,ratio_imag,rho,quasi_energy,projected_residual")
            for k in ks,c in s.cases
                modes=checked_spectrum_call(()->Base.invokelatest(sector_modes,s.op,s.r,c.right,c.sector,k;num,tol=stol,deadline))
                for (j,(v,res)) in enumerate(zip(modes.values,modes.residuals))
                    rho=abs(v);println(io,join((c.name,k,j,real(v),imag(v),rho,-log(rho),res),','))
                    count+=1;maxres=max(maxres,res)
                end
                flush(io);println("SPECTRUM panel=$panel k/pi=$(round(k/pi;digits=3)) sector=$(c.name) dim=$(modes.dimension) elapsed=$(round(time()-START;digits=1))");flush(stdout)
            end
        end
    catch err
        write_json(joinpath(out,"summary.json"),(;status="error",count,max_projected_residual=maxres,error=sprint(showerror,err),wall_seconds=time()-START))
        rethrow()
    end
    write_json(joinpath(out,"summary.json"),(;status="finished",count,max_projected_residual=maxres,wall_seconds=time()-START,peak_rss_bytes=Sys.maxrss()))
end
end
if abspath(PROGRAM_FILE)==@__FILE__;FilteredFigure7.main();end
