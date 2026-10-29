struct BoundaryPair
    r0::InfiniteMPS
    r1::InfiniteMPS
    env0
    env1
    adjoint0::InfiniteMPS
    adjoint1::InfiniteMPS
    adjoint_env0
    adjoint_env1
    diagnostics::NamedTuple
end
function product_boundary(net;phase=1)
    W=space(net.ket_action,1); b=Z2Space(0=>1)
    d=dim(W);D=isqrt(d)
    f=reshape(convert(Array,net.fusion),d,d)
    a=phase*normalize(f*vec(Matrix{Float64}(I,D,D)))
    InfiniteMPS([TensorMap(reshape(ComplexF64.(a),1,d,1),b⊗W←b)])
end
function apply_group(psi,g)
    tensors=map(1:length(psi)) do i
        a=psi.AL[i]
        @tensor t[l p; r] := a[l q; r]*g[p; q]
        t
    end
    InfiniteMPS(tensors)
end
boundary_overlap(a,b)=Float64(abs(dot(a,b))/sqrt(abs(dot(a,a)*dot(b,b))))
distinct_boundaries(a,b;tol=1e-8)=boundary_overlap(a,b)<1-tol
remove_zero_support(a)=changebonds(a,SvdCut(;trscheme=trunctol(;atol=1e-14)))
function boundary_guess(net,chi)
    chi==1 && return product_boundary(net)
    W=space(net.ket_action,1); b=Z2Space(0=>chi÷2,1=>chi÷2)
    rng=MersenneTwister(2026)
    a=randn(rng,ComplexF64,b⊗W←b)
    InfiniteMPS([a])
end
function adjoint_transfer(op)
    tensors=map(1:length(op)) do i
        a=op[i]; W=space(a,1)
        TensorMap(permutedims(conj(convert(Array,a)),(1,3,2,4)),W⊗W←W⊗W;tol=1e-12)
    end
    InfiniteMPO(tensors)
end
density(psi,op,env)=ComplexF64(prod(expectation_value(psi,op,env)))
boundary_solve(guess,op;kwargs...)=checked_spectrum_call(()->_boundary_solve(guess,op;kwargs...))
function _boundary_solve(guess,op;tol,maxiter,deadline)
    finalize=(iter,psi,h,env)->begin
        time()<deadline || error("boundary time budget reached")
        if iter%10==0
            println("BOUNDARY iter=$iter projected_residual=$(MPSKit.calc_galerkin(psi,h,psi,env))");flush(stdout)
        end
        psi,env
    end
    # Leave room for the support truncation and rebuilt environment below;
    # their residual is the one tested against the requested tolerance.
    alg=VUMPS(;tol=tol/10,maxiter,verbosity=0,alg_eigsolve=MPSKit.Defaults.alg_eigsolve(;ishermitian=false),finalize)
    psi,_,_=leading_boundary(guess,op,alg)
    psi=remove_zero_support(psi)
    env=environments(psi,op;tol)
    residual=Float64(MPSKit.calc_galerkin(psi,op,psi,env))
    psi,env,residual
end
function boundary_chis(options)
    chis=Tuple(get(options,:chis,(1,2,4,8,16,32)))
    !isempty(chis) && all(c->c isa Integer && (c==1 || c>0 && iseven(c)),chis) ||
        throw(ArgumentError("boundary dimensions must be 1 or positive even integers"))
    return chis
end
solve_boundaries(net;options=(;))=checked_spectrum_call(()->_solve_boundaries(net;options))
function boundary_converged(d;tol=1e-10,fixedpoint_tol=tol)
    all(x->isfinite(x) && x>0,(tol,fixedpoint_tol)) || throw(ArgumentError("boundary tolerances must be finite and positive"))
    local_error=maximum((d.projected_residual0,d.projected_residual1,d.adjoint_residual0,d.adjoint_residual1,
        d.density_mismatch,1-d.diagonal_overlap0,1-d.diagonal_overlap1))
    approximation_error=max(d.adjoint_density_mismatch,d.fidelity_defect)
    local_error<=tol && approximation_error<=fixedpoint_tol && d.overlap01<1-1e-8 && d.adjoint_right_overlap>1e-8
end
function _solve_boundaries(net;options=(;))
    tol=get(options,:tol,1e-10); maxiter=get(options,:maxiter,100); deadline=get(options,:deadline,Inf)
    fixedpoint_tol=get(options,:fixedpoint_tol,tol)
    last=nothing
    for chi in boundary_chis(options)
        r0,env0,res0=boundary_solve(boundary_guess(net,chi),net.mpo;tol,maxiter,deadline)
        r1=apply_group(r0,net.ket_action); env1=environments(r1,net.mpo;tol)
        res1=Float64(MPSKit.calc_galerkin(r1,net.mpo,r1,env1))
        lambda0=density(r0,net.mpo,env0);lambda1=density(r1,net.mpo,env1)
        mismatch=abs(lambda0-lambda1)/abs(lambda0)
        overlap01=boundary_overlap(r0,r1)
        C=net.ket_action*net.bra_action
        c0=boundary_overlap(r0,apply_group(r0,C));c1=boundary_overlap(r1,apply_group(r1,C))
        opadj=adjoint_transfer(net.mpo)
        a0,aenv0,ares0=boundary_solve(boundary_guess(net,chi),opadj;tol,maxiter,deadline)
        a1=apply_group(a0,net.ket_action); aenv1=environments(a1,opadj;tol)
        ares1=Float64(MPSKit.calc_galerkin(a1,opadj,a1,aenv1))
        alambda=density(a0,opadj,aenv0); amismatch=abs(conj(alambda)-lambda0)/abs(lambda0)
        # Fidelity per site is a separate fixed-point diagnostic, not a norm bound.
        fidelity_defect=max(0.0,1-boundary_overlap(r0,net.mpo*r0))
        ar_overlap=boundary_overlap(a0,r0)
        diagnostics=(;chi=dim(space(r0.AL[1],1)),projected_residual0=res0,projected_residual1=res1,
            adjoint_residual0=ares0,adjoint_residual1=ares1,density0_real=real(lambda0),density0_imag=imag(lambda0),
            density1_real=real(lambda1),density1_imag=imag(lambda1),density_mismatch=mismatch,
            adjoint_density_mismatch=amismatch,overlap01,diagonal_overlap0=c0,diagonal_overlap1=c1,
            single_layer_overlap=boundary_overlap(r1,apply_group(r0,net.ket_action)),fidelity_defect,adjoint_right_overlap=ar_overlap,
            solver_tol=tol,fixedpoint_tol)
        converged=boundary_converged(diagnostics;tol,fixedpoint_tol)
        diagnostics=merge(diagnostics,(;converged))
        last=BoundaryPair(r0,r1,env0,env1,a0,a1,aenv0,aenv1,diagnostics)
        converged && return last
    end
    error("no reliable distinct boundary pair at requested dimensions: $(last.diagnostics)")
end
