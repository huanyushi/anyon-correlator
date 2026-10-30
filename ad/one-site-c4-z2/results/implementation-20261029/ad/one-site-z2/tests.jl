using Test, LinearAlgebra, Random
using TensorKit
include("core.jl")
using .OneSiteZ2AD
const A=OneSiteZ2AD
function torus_wavefunction(t)
    @tensor state[p,q,r,s] := t[p,i,j,k,l]*t[q,m,l,n,j]*t[r,k,o,i,u]*t[s,n,u,m,o]
    return vec(state)/norm(state)
end
@testset "numerical failures have typed checkpoint-safe outcomes" begin
    b=A.ParameterBasis(2);x=A.initial_parameters(b,9)
    @test_throws A.Legacy.GradientFailure A.parameter_pullback(x,b,fill(NaN,length(A.split_parameters(x,b))))
    @test_throws A.Legacy.GradientFailure A.checked_legacy_gradient(()->error("AD parameter gradient is nonfinite"))
    @test_throws ErrorException A.checked_legacy_gradient(()->error("unrelated programming defect"))
    bad=copy(x);bad[end]=10000
    @test_throws A.Legacy.ContractionFailure A.checked_trial(bad,b)
    mktempdir() do out
        result=A.optimize_steps(x,b;outdir=joinpath(out,"run"),
            evaluate=(z;gradient=true)->A.checked_legacy_gradient(()->error("AD parameter gradient is nonfinite")))
        @test result.status=="gradient_failed"
        saved,_,_=A.load_state(joinpath(out,"run","final.json"))
        @test saved==x
    end
end
@testset "audit outcomes finalize authoritative state and summary" begin
    b=A.ParameterBasis(2);x=A.initial_parameters(b,6)
    cases=[("stationary",true), ("gradient_unstable",false), ("audit_rejected",false)]
    mktempdir() do dir
        for (status,accepted) in cases
            out=joinpath(dir,status);mkpath(out)
            A.save_state(joinpath(out,"final.json"),x,b;metadata=(;accepted=false))
            result=A.finalize_run((;status="stationary_candidate",steps=0,accepted=false),x,b,Dict(),out;
                deadline=time()+30,assess=()->(;status,accepted,reasons=accepted ? [] : [status]))
            @test result.accepted==accepted && result.status==status
            @test A.JSON.parsefile(joinpath(out,"summary.json"))["accepted"]==accepted
            @test last(A.load_state(joinpath(out,"final.json")))["accepted"]==accepted
        end
        for (name,deadline,assess) in (("expired",time()-1,()->error("must not run")),
            ("failed",time()+30,()->throw(A.Legacy.GradientFailure("failed audit"))))
            out=joinpath(dir,name);mkpath(out)
            A.save_state(joinpath(out,"final.json"),x,b;metadata=(;accepted=false))
            result=A.finalize_run((;status="stationary_candidate",steps=0,accepted=false),x,b,Dict(),out;deadline,assess)
            @test !result.accepted
            @test result.status==(name=="expired" ? "budget" : "audit_failed")
        end
    end
end
function north_east_block(st)
    D=length(st.lambda);out=zeros(2,2,D,D,D,D)
    for p in 1:2,q in 1:2,n in 1:D,e in 1:D,s in 1:D,w in 1:D,a in 1:D,b in 1:D
        out[p,q,n,e,s,w]+=st.V[b,a,s,w]*st.P[p,a,e]*st.P[q,b,n]*
            st.lambda[a]*st.lambda[b]*st.lambda[e]*st.lambda[n]
    end
    reshape(out,4,D,D,D,D)
end
@testset "closed-network gauge equivalence and exact stabilizers" begin
    for D in (2,4)
        b=A.ParameterBasis(D);x=A.initial_parameters(b,29)
        x[end-D+1:end].=range(-0.4,0.6;length=D)
        balanced=A.Legacy.blocked_arrays(A.split_parameters(x,b),D)[1,1]
        @test torus_wavefunction(balanced)≈torus_wavefunction(north_east_block(A.tensors(x,b))) atol=1e-12
        exact=A.initial_parameters(b,1;exact=true)
        psi=torus_wavefunction(A.Legacy.blocked_arrays(A.split_parameters(exact,b),D)[1,1])
        for sites in ((1,3,2,6),(3,1,4,8),(5,7,6,2),(7,5,8,4))
            mask=sum(1<<(i-1) for i in sites)
            @test psi≈psi[[xor(i,mask)+1 for i in 0:255]] atol=1e-12
        end
        for sites in ((1,5,2,4),(3,7,4,2),(5,1,6,8),(7,3,8,6))
            signs=[(-1)^sum((i>>(s-1))&1 for s in sites) for i in 0:255]
            @test psi≈psi.*signs atol=1e-12
        end
    end
end
@testset "charge and C4 parameterization at D=2,4" begin
    for D in (2,4)
        basis=A.ParameterBasis(D)
        x=A.initial_parameters(basis,19)
        st=A.tensors(x,basis)
        @test length(x)==(D==2 ? 10 : 54)
        @test st.V≈permutedims(st.V,(2,3,4,1)) atol=1e-14
        @test st.P≈permutedims(st.P,(1,3,2)) atol=1e-14
        @test all(st.V[i]==0 for i in CartesianIndices(st.V) if isodd(sum(basis.charges[j] for j in Tuple(i))))
        @test all(st.P[p,a,b]==0 for p in 1:2,a in 1:D,b in 1:D if basis.charges[a]!=basis.charges[b])
        @test x==A.initial_parameters(basis,19)
        @test all(>(0),st.lambda)
        @test_throws ArgumentError A.tensors(zeros(length(x)),basis)
        y=A.split_parameters(x,basis)
        @test length(y)==A.Legacy.parameter_count(D;uniform=true)
        old=A.Legacy.split_parameters(y,D)
        @test all(v->v≈st.V,old.V)
        @test all(p->p≈st.P,[old.E;old.N])
        @test all(w->w≈st.lambda,old.weights)
        coeff=randn(MersenneTwister(4),length(y))
        f(z)=dot(A.split_parameters(z,basis),coeff)
        g=only(A.Zygote.gradient(f,x));direction=normalize(randn(MersenneTwister(5),length(x)))
        @test dot(g,direction)≈(f(x+1e-5direction)-f(x-1e-5direction))/2e-5 atol=1e-8 rtol=1e-6
    end
end
@testset "portable AD provenance and exact embedding" begin
    for D in (2,4)
        b=A.ParameterBasis(D);x=A.initial_parameters(b,1;exact=true)
        st=A.tensors(x,b)
        active=[1,findfirst(==(1),b.charges)]
        @test st.V[active,active,active,active]≈A.tensors(A.initial_parameters(A.ParameterBasis(2),1;exact=true),A.ParameterBasis(2)).V
        mktempdir() do dir
            path=joinpath(dir,"state.json")
            A.save_state(path,x,b;metadata=(;origin="exact_zero_initialization",accepted=false))
            y,c,md=A.load_state(path)
            @test y==x && c.charges==b.charges
            @test md["accepted"]===false
            @test_throws ArgumentError A.save_state(path,x,b;metadata=(;))
            data=A.JSON.parsefile(path);data["algorithm"]="simple_update"
            write(path,A.JSON.json(data))
            @test_throws ArgumentError A.load_state(path)
        end
    end
end
@testset "shared-parameter chain rule" begin
    b=A.ParameterBasis(4);x=A.initial_parameters(b,11)
    y=A.split_parameters(x,b);g=randn(MersenneTwister(2),length(y))
    actual=A.parameter_pullback(x,b,g)
    d=normalize(randn(MersenneTwister(3),length(x)))
    @test dot(actual,d)≈dot(g,(A.split_parameters(x+1e-5d,b)-A.split_parameters(x-1e-5d,b))/2e-5) atol=1e-8 rtol=1e-6
end
@testset "budget preserves the last complete state" begin
    b=A.ParameterBasis(2);x=A.initial_parameters(b,4)
    # A real scalar quadratic isolates optimizer control flow from CTMRG cost.
    evaluate(z;gradient=true)= (;energy=dot(z,z)/2,gradient=gradient ? copy(z) : nothing)
    mktempdir() do dir
        result=A.optimize_steps(x,b;evaluate,outdir=joinpath(dir,"expired"),deadline=time()-1)
        @test result.status=="budget" && result.steps==0 && !result.accepted
        saved,_,_=A.load_state(joinpath(dir,"expired","final.json"))
        @test saved==x
        result=A.optimize_steps(x,b;evaluate,outdir=joinpath(dir,"steps"),maxsteps=2)
        @test result.status=="step_limit" && result.steps==2 && !result.accepted
        saved,_,_=A.load_state(joinpath(dir,"steps","final.json"))
        @test dot(saved,saved)<dot(x,x)
        @test_throws ArgumentError A.optimize_steps(x,b;evaluate,outdir=joinpath(dir,"steps"))
        failed(z;gradient=true)=throw(A.Legacy.GradientFailure("test failed adjoint"))
        result=A.optimize_steps(x,b;evaluate=failed,outdir=joinpath(dir,"failed"))
        @test result.status=="gradient_failed" && !result.accepted
    end
end
@testset "production gradient evidence and continuation identity" begin
    b=A.ParameterBasis(4)
    p=Dict("passed"=>true,"D"=>4,"charges"=>b.charges,"chi"=>8,"hx"=>0.,"hz"=>0.1,
        "parameter_count"=>54,"gradient_method"=>"unrolled","unroll_steps"=>20)
    @test A.probe_matches(p,b,0.,0.1;chi=8)
    @test !A.probe_matches(p,b,0.,0.2;chi=8)
    @test !A.probe_matches(p,b,0.,0.1;chi=16)
    @test !A.probe_matches(p,A.ParameterBasis(4;n_even=1),0.,0.1;chi=8)
    @test_throws ArgumentError A.physical_evaluator(b,0.,0.1;chi=16,probe=p)
    md=Dict("hx"=>0.,"hz"=>0.1,"origin"=>"symmetric_random","seed"=>19)
    @test A.check_resume(md,0.,0.1)==md
    @test_throws ArgumentError A.check_resume(md,0.,0.2)
    mktempdir() do dir
        probe=joinpath(dir,"probe.json");write(probe,A.JSON.json(p))
        out=joinpath(dir,"run")
        summary=A.run(Dict("D"=>"4","probe"=>probe,"out"=>out,"budget"=>"0","init"=>"random"))
        @test summary.status=="budget"
        _,_,saved=A.load_state(joinpath(out,"final.json"))
        @test saved["origin"]=="symmetric_random" && saved["hx"]==0. && saved["hz"]==0.1
        @test_throws ArgumentError A.run(Dict("resume"=>joinpath(out,"final.json"),"hz"=>"0.2","out"=>joinpath(dir,"bad"),"probe"=>probe))
    end
end
