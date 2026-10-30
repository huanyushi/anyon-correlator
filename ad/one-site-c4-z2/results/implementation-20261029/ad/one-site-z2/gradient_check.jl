# Deliberately separate from the lightweight software tests.
using LinearAlgebra, Random
include("core.jl")
const A=OneSiteZ2AD
D=parse(Int,get(ENV,"AD_D","2"));chi=parse(Int,get(ENV,"AD_CHI","8"))
b=A.ParameterBasis(D)
x=A.initial_parameters(b,71;exact=true)+0.03A.initial_parameters(b,71)
out=get(ENV,"AD_OUT","")
if !isempty(out)
    ispath(out) && error("output directory already exists: $out")
    mkpath(out)
    A.save_state(joinpath(out,"initial.json"),x,b;metadata=(;
        origin="perturbed_exact_zero_initialization",seed=71,amplitude=0.03,
        hx=0.,hz=0.1,chi,accepted=false,purpose="gradient_validation"))
end
H,_=A.Legacy.operators(0.,0.1)
println("contract D=$D chi=$chi parameters=$(length(x))");flush(stdout)
r=A.contract(x,b,H;chi,tol=1e-10,maxiter=500)
println("energy=$(r.energy) residual=$(r.info.convergence_error)");flush(stdout)
gradients=[];rows=[]
for n in (20,40)
    t=time();value,g=A.energy_gradient(x,b,H,r;method=:unrolled,unroll_steps=n)
    push!(gradients,g)
    println("N=$n energy=$value norm=$(norm(g)) seconds=$(time()-t)");flush(stdout)
    d=normalize(g)
    for h in (1e-4,1e-5,1e-6)
        fd=(A.unrolled_energy(x+h*d,b,H,r;steps=n)-A.unrolled_energy(x-h*d,b,H,r;steps=n))/(2h)
        err=abs(fd-dot(g,d));pass=err<=max(1e-7,1e-3abs(dot(g,d))) && abs(fd)>1e-7
        push!(rows,(;n,h,fd,ad=dot(g,d),err,pass))
        if !isempty(out)
            open(joinpath(out,"directions.jsonl"),"a") do io
                println(io,A.JSON.json(last(rows)));flush(io)
            end
        end
        println(last(rows));flush(stdout)
    end
end
delta=norm(gradients[1]-gradients[2])
passed=all(row.pass for row in rows) && delta<=2e-6
println("gradient_delta=$delta passed=$passed");flush(stdout)
if !isempty(out)
    open(joinpath(out,"summary.json"),"w") do io
        A.JSON.print(io,(;D,charges=b.charges,chi,hx=0.,hz=0.1,parameter_count=length(x),
            gradient_method="unrolled",unroll_steps=20,delta,passed,rows));flush(io)
    end
end
passed || error("gradient validation failed; optimization is not authorized")
