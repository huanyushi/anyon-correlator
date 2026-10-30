module FilteredToricBenchmark
const START=time()
include("../../transfer-spectrum/numerics.jl")
using JSON
include("tensor.jl")

write_json(path,data)=open(io->(JSON.print(io,data);println(io)),path,"w")

function options(args)
    opts=Dict("out"=>"","chi"=>"8","nk"=>"41","num"=>"24","tol"=>"1e-6",
        "spectrum-tol"=>"1e-7","budget"=>"1800")
    for arg in args
        k,v=split(replace(arg,r"^--"=>""),'=';limit=2)
        haskey(opts,k) || error("unknown option $k")
        opts[k]=v
    end
    isempty(opts["out"]) && error("--out=NEW_DIRECTORY is required")
    opts
end

function main(args=ARGS)
    opts=options(args);out=opts["out"]
    chi=parse(Int,opts["chi"]);nk=parse(Int,opts["nk"]);num=parse(Int,opts["num"])
    tol=parse(Float64,opts["tol"]);spectrum_tol=parse(Float64,opts["spectrum-tol"])
    budget=parse(Float64,opts["budget"]);deadline=START+budget
    boundary_chis((;chis=(chi,)))
    nk>=2 && num>0 && all(x->isfinite(x) && x>0,(tol,spectrum_tol,budget)) || error("invalid numerical options")
    ispath(out) && error("output already exists: $out")
    mkpath(out)
    # Save only this driver's actual sources and the fixed package environment.
    root=normpath(joinpath(@__DIR__,"../.."))
    sources=vcat(["transfer-spectrum/$f.jl" for f in ("numerics","network","boundary","spectrum")],
        ["benchmarks/filtered-toric-code/$f" for f in ("run.jl","tensor.jl","plot.py")],
        ["julia-env/Project.toml","julia-env/Manifest.toml"])
    for file in sources
        dest=joinpath(out,"source",file);mkpath(dirname(dest));cp(joinpath(root,file),dest)
    end
    ks=collect(range(-pi,pi;length=nk))
    write_json(joinpath(out,"settings.json"),(;opts,betax=.45,betaz=.55,chi,nk,num,tol,spectrum_tol,
        momenta=ks,geometry="four-spin diagonal",horizontal_vector=[1,1],vertical_vector=[-1,1],
        domain_wall_background_phase=pi,julia=string(VERSION),mpskit=string(pkgversion(MPSKit))))
    count=0;maxres=0.;status="finished"
    try
        # 1. Analytic PEPS -> double-layer ordinary row MPO.
        A=filtered_tensor();net=build_transfer_from_tensor(A,[0,1])
        write_json(joinpath(out,"tensor.json"),(;shape=size(A),A=vec(A),charges=[0,1]))

        # 2. One VUMPS solve. The other background follows from the virtual Z2.
        println("BOUNDARY chi=$chi tol=$tol");flush(stdout)
        r0,env,res=boundary_solve(boundary_guess(net,chi),net.mpo;tol,maxiter=300,deadline)
        r1=apply_group(r0,net.ket_action)
        # This sign gives the paper's convention: the flux minimum is at k=π.
        r1=InfiniteMPS([-r1.AL[1]],[-r1.AR[1]],[r1.C[1]],[-r1.AC[1]])
        a=convert(Array,r0.AL[1])
        write_json(joinpath(out,"boundary.json"),(;chi,shape=size(a),AL_real=vec(real.(a)),AL_imag=vec(imag.(a)),
            charge_dims=[dim(space(r0.AL[1],1),Z2Irrep(q)) for q in 0:1],residual=res,
            density_real=real(density(r0,net.mpo,env)),density_imag=imag(density(r0,net.mpo,env))))

        # 3. At each momentum, diagonalize the same QP effective operator in
        #    four sectors: same/different background and even/odd charge.
        open(joinpath(out,"spectrum.csv"),"w") do io
            println(io,"g,q,k,branch,ratio_real,ratio_imag,rho,quasi_energy,projected_residual")
            flush(io)
            for k in ks,g in 0:1,q in 0:1
                time()<deadline || error("time budget reached")
                right=g==0 ? r0 : r1
                modes=solve_modes(net.mpo,r0,right,q,k;num,options=(;tol=spectrum_tol,deadline))
                for (branch,(r,residual)) in enumerate(zip(modes.ratios,modes.residuals))
                    rho=abs(r);epsilon=rho>0 ? -log(rho) : Inf
                    println(io,join((g,q,k,branch,real(r),imag(r),rho,epsilon,residual),','))
                    count+=1;maxres=max(maxres,residual)
                end
                flush(io)
                println("SPECTRUM k/pi=$(round(k/pi;digits=3)) g=$g q=$q modes=$(length(modes.ratios)) elapsed=$(round(time()-START;digits=1))");flush(stdout)
            end
        end
    catch err
        status=time()>=deadline ? "budget" : "error"
        write_json(joinpath(out,"summary.json"),(;status,count,max_projected_residual=maxres,error=sprint(showerror,err),wall_seconds=time()-START))
        rethrow()
    end
    write_json(joinpath(out,"summary.json"),(;status,count,max_projected_residual=maxres,wall_seconds=time()-START,peak_rss_bytes=Sys.maxrss()))
end
end
if abspath(PROGRAM_FILE)==@__FILE__;FilteredToricBenchmark.main();end
