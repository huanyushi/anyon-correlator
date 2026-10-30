using Test, LinearAlgebra, TensorKit
BLAS.set_num_threads(1)
include("core.jl")
using .OneSiteZ2SU
const S=OneSiteZ2SU

group=isempty(ARGS) ? "all" : only(ARGS)
group in ("all","--group=tensors","--group=updates","--group=runner") || error("unknown group $group")

if group in ("all","--group=tensors")
    @testset "Z₂ representation and spatial constraints" begin
        st=random_state(5); a=dense_arrays(exact_reference_state())
        @test dim(virtual_space(),Z2Irrep(0))==dim(virtual_space(),Z2Irrep(1))==1
        @test dim(physical_space(),Z2Irrep(0))==2 && dim(physical_space(),Z2Irrep(1))==0
        @test maximum(values(symmetry_residuals(st)))<=1e-12
        @test dense_arrays(random_state(5))==dense_arrays(st)
        @test norm(dense_arrays(st).P-a.P)>1e-3
        @test a.V ≈ [iseven(n+e+s+w) ? 1/sqrt(8) : 0 for n in 0:1,e in 0:1,s in 0:1,w in 0:1] atol=1e-12
        @test a.P ≈ [u==v ? (-1.0)^(p*u)/2 : 0 for p in 0:1,u in 0:1,v in 0:1] atol=1e-12
        @test a.lambda ≈ [1,1]/sqrt(2)
        v=zeros(2,2,2,2); v[1,1,2,2]=1
        projected=state_from_arrays(v,a.P,a.lambda)
        @test maximum(values(symmetry_residuals(projected)))<=1e-12
        @test norm(dense_arrays(projected).V-permutedims(dense_arrays(projected).V,(1,3,2,4)))>0.1
    end
end

if group in ("all","--group=updates")
    include("../split/core.jl")
    @testset "Local contractions and charge-preserving updates" begin
        st=random_state(5); a=dense_arrays(st); w=SplitSU.initial_state(:exact,0)
        for n in keys(w.tensors); w.tensors[n]=copy(n[1]==:V ? a.V : a.P); end
        for b in keys(w.weights); w.weights[b]=copy(a.lambda); end
        @test reshape(convert(Array,S.star_tensor(st)),4,4,4,4) ≈ first(SplitSU.star_cluster(w,1,1)) atol=1e-12
        @test vec(convert(Array,S.plaquette_tensor(st))) ≈ vec(SplitSU.ring_tensor(first(SplitSU.ring_cluster(w,1,1)))) atol=1e-12
        before=dense_arrays(st)
        @test_throws ArgumentError sweep!(st,0.05;als_sweeps=0)
        @test dense_arrays(st)==before
        d=sweep!(st,0.05;als_sweeps=12)
        @test isfinite(d.star.discarded_weight) && isfinite(d.plaquette.fit_residual)
        @test max(d.star.gauge_residual,d.plaquette.gauge_residual)<=1e-8
        @test maximum(values(symmetry_residuals(st)))<=1e-12
        @test all(isfinite,values(local_residuals(st)))
        ref=exact_reference_state(); q=virtual_space(); z=TensorMap([1.0 0;0 -1],q←q)
        ps=[S.outgoing(ref.P) for _ in 1:4]; ps[1]=S.mode(ps[1],z,2)
        vv,pp=S.align_star(ref,S.mode(ref.V,z,1),ps)
        @test norm(S.canonical_v(vv)-ref.V)<1e-12
        @test all(p->norm(S.canonical_p(p)-ref.P)<1e-12,pp)
        a=dense_arrays(ref); zero_support=state_from_arrays(a.V,a.P,[1.0,0])
        @test all(isfinite,convert(Array,S.support_inverse(zero_support.lambda)))
        @test maximum(values(local_residuals(ref)))<1e-12
        sweep!(ref,0.05;als_sweeps=12)
        @test maximum(values(local_residuals(ref)))<1e-12
    end
    @testset "X-polarized ALS stationary point is not a ground state" begin
        a=dense_arrays(random_state(5))
        st=state_from_arrays(a.V,[u==v ? 0.5 : 0.0 for p in 0:1,u in 0:1,v in 0:1],a.lambda)
        d=plaquette_step!(st,0.05;als_sweeps=12)
        expected=sinh(0.05)^2/(cosh(0.05)^2+sinh(0.05)^2)
        @test d.fit_residual ≈ expected atol=1e-12
        @test local_residuals(st).plaquette ≈ sqrt(2) atol=1e-12
        good=(;energy=-1.,vals=ones(8),converged=true,residual=0.,chi=8)
        @test !S.audit([good,merge(good,(chi=16,))],st,"finished").passed
    end
end

if group in ("all","--group=runner")
    @testset "State I/O and same-trajectory continuation" begin
        mktempdir() do dir
            st=random_state(5); file=joinpath(dir,"state.json")
            md=(;origin="numerical_random_optimization",seed=5,step=300,dt=0.05,als_sweeps=12,accepted=false)
            save_state(file,st;metadata=md)
            resumed,step=S.resume_state(file;seed=5,dt=0.05,als_sweeps=12)
            @test dense_arrays(resumed)==dense_arrays(st) && step==300
            @test_throws ArgumentError S.resume_state(file;seed=6,dt=0.05,als_sweeps=12)
            @test_throws ArgumentError S.resume_state(file;seed=5,dt=0.1,als_sweeps=12)
            @test_throws ArgumentError save_state(file,st;metadata=md)
            @test_throws ArgumentError S.prepare_output(dir)
            bad=joinpath(dir,"bad.json"); data=S.JSON.parsefile(file); data["V"][9]=1e-9
            S.write_json(bad,data)
            @test_throws ArgumentError load_state(bad)
            ref=joinpath(dir,"ref.json"); save_state(ref,exact_reference_state();metadata=merge(md,(origin="reference",)))
            @test_throws ArgumentError S.resume_state(ref;seed=5,dt=0.05,als_sweeps=12)
        end
    end
    @testset "Physical acceptance and hard time budget" begin
        st=exact_reference_state(); good=(;energy=-1.,vals=ones(8),converged=true,residual=0.,chi=8)
        other=merge(good,(chi=16,))
        @test S.audit([good,other],st,"finished").passed
        @test !S.audit([good,other],st,"budget").passed
        for bad in ((converged=false,),(vals=zeros(8),),(energy=NaN,))
            @test !S.audit([merge(good,bad),other],st,"finished").passed
        end
        mktempdir() do dir
            marker=joinpath(dir,"started")
            @test S.supervise(`touch $marker`,dir;budget=2.,deadline=time()-1)=="budget"
            @test !isfile(marker)
            save_state(joinpath(dir,"latest.json"),random_state(5);metadata=Dict())
            @test S.supervise(`sleep 5`,dir;budget=0.05)=="budget"
            @test S.JSON.parsefile(joinpath(dir,"summary.json"))["passed"]==false
            @test dense_arrays(first(load_state(joinpath(dir,"latest.json"))))==dense_arrays(random_state(5))
            open(joinpath(dir,"worker.log"),"w") do io
                @test S.supervise(`printf progress`,dir;budget=2.,output=io)=="finished"
            end
            @test read(joinpath(dir,"worker.log"),String)=="progress"
        end
    end
end
