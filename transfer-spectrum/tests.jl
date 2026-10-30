# Compact implementation checks, separate from physical spectrum production.
using Test, LinearAlgebra, TensorKit
include("core.jl")
using .ToricTransferSpectrum
const T=ToricTransferSpectrum
const S=T.OneSiteZ2SU
const A=T.OneSiteZ2AD

@testset "V/P/lambda input and virtual symmetry" begin
    st=S.state_from_arrays(S.dense_arrays(S.random_state(5)).V,
        S.dense_arrays(S.random_state(5)).P,[0.8,0.6])
    a=S.dense_arrays(st);net=build_transfer(st);expected=zeros(4,2,2,2,2)
    for p in 0:1,q in 0:1,n in 1:2,e in 1:2,s in 1:2,w in 1:2,b in 1:2,c in 1:2
        expected[p+2q+1,n,e,s,w]+=a.V[c,b,s,w]*a.P[p+1,b,e]*a.P[q+1,c,n]*a.lambda[b]*a.lambda[e]*a.lambda[c]*a.lambda[n]
    end
    @test net.A≈expected atol=1e-12
    @test net.step==1
    @test length(net.mpo)==1
    for g in (net.ket_action,net.bra_action)
        @test g*g≈id(space(g,1)) atol=1e-12
        q=copy(net.Q)
        for i in 1:4;q=T.mode_array(q,convert(Array,g),i);end
        @test norm(q-net.Q)/norm(net.Q)<1e-12
    end
    b=A.ParameterBasis(4)
    net4=build_transfer(merge(A.tensors(A.initial_parameters(b,1;exact=true),b),(;charges=b.charges)))
    @test size(net4.A)==(4,4,4,4,4)
    @test size(net4.Q)==(16,16,16,16)
end

@testset "sector labels, residuals and frozen input provenance" begin
    @test [SectorLabel(g,q).name for g in 0:1 for q in 0:1]==["1","e","m","epsilon"]
    @test_throws ArgumentError SectorLabel(2,0)
    @test T.mode_status(0.5,1e-9,4;tol=1e-7)=="converged"
    @test T.mode_status(0.5,1e-3,4;tol=1e-7)=="not_converged"
    @test T.mode_status(1e-15,1e-10,4;tol=1e-7)=="below_resolution"
    @test T.boundary_chis((;chis=(8,12)))==(8,12)
    @test_throws ArgumentError T.boundary_chis((;chis=(3,)))
    mktempdir() do dir
        b=A.ParameterBasis(4);x=A.initial_parameters(b,8)
        path=joinpath(dir,"candidate.json")
        A.save_state(path,x,b;metadata=(;accepted=false,status="stationary_candidate",
            origin="perturbed_exact_zero_initialization",seed=8,hx=0.,hz=0.1))
        @test_throws ArgumentError T.load_production_state(path)
        st,md=T.load_production_state(path;reference=true)
        @test st.V≈A.tensors(x,b).V
        @test !md["accepted"]
    end
end
