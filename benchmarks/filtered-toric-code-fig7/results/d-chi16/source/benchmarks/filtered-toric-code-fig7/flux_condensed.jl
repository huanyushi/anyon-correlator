# Both layer charges are conserved independently in the unique symmetric fixed point.
function flux_condensed_setup(A,chi;tol,deadline)
    chi%4==0 || error("flux-condensed chi must be divisible by 4")
    net=build_transfer_from_tensor(A,[0,1])
    S=ProductSector{Tuple{Z2Irrep,Z2Irrep}}
    W=Vect[S]((0,0)=>1,(0,1)=>1,(1,0)=>1,(1,1)=>1)
    b=Vect[S]((0,0)=>chi÷4,(0,1)=>chi÷4,(1,0)=>chi÷4,(1,1)=>chi÷4)
    # Raw double-layer order is ket-fast; product irreps are lexicographic.
    Q=net.Q_raw[[1,3,2,4],[1,3,2,4],[1,3,2,4],[1,3,2,4]]/net.normalization
    op=InfiniteMPO([TensorMap(permutedims(Q,(4,3,1,2)),W⊗W←W⊗W;tol=1e-12)])
    zk=TensorMap(Matrix(Diagonal([1.,1,-1,-1])),W←W)
    zb=TensorMap(Matrix(Diagonal([1.,-1,1,-1])),W←W)
    guess=InfiniteMPS([randn(MersenneTwister(2026),ComplexF64,b⊗W←b)])
    r,env,res=boundary_solve(guess,op;tol,maxiter=400,deadline)
    cases=[(;name="q$(qk)$(qb)",right=r,sector=S((qk,qb)),color="blue",
        marker=qk!=qb ? "circle" : qk==0 ? "plus" : "square") for qk in 0:1 for qb in 0:1]
    (;op,r,env,res,cases,actions=[zk,zb],basis="Z2-times-Z2")
end
