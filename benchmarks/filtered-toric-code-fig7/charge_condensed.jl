# No imposed boundary symmetry: select one of the four broken-symmetry fixed points.
function charge_condensed_setup(A,chi;tol,deadline)
    net=build_transfer_from_tensor(A,[0,1]);W=ComplexSpace(4);b=ComplexSpace(chi)
    op=InfiniteMPO([TensorMap(ComplexF64.(permutedims(net.Q_raw/net.normalization,(4,3,1,2))),W⊗W←W⊗W)])
    zk=TensorMap(Matrix{ComplexF64}(Diagonal([1,-1,1,-1])),W←W)
    zb=TensorMap(Matrix{ComplexF64}(Diagonal([1,1,-1,-1])),W←W)
    rng=MersenneTwister(2026)
    a=.02randn(rng,ComplexF64,chi,4,chi)
    a[1,:,1].+=.5 # X-positive ket and bra selects a symmetry-broken boundary.
    guess=InfiniteMPS([TensorMap(a,b⊗W←b)])
    r,env,res=boundary_solve(guess,op;tol,maxiter=400,deadline)
    cases=[(;name="identity",right=r,sector=Trivial(),color="blue",marker="dot"),
        (;name="ket-string",right=phase_flip(apply_group(r,zk)),sector=Trivial(),color="red",marker="dot"),
        (;name="ket-bra-string",right=apply_group(r,zk*zb),sector=Trivial(),color="green",marker="dot")]
    (;op,r,env,res,cases,actions=[zk,zb],basis="dense-ket-bra")
end
