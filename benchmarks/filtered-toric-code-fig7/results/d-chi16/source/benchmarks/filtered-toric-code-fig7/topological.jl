# The diagonal Z2 charge is conserved; a single-layer action changes the fixed point.
function topological_setup(A,chi;tol,deadline)
    net=build_transfer_from_tensor(A,[0,1])
    r,env,res=boundary_solve(boundary_guess(net,chi),net.mpo;tol,maxiter=400,deadline)
    other=apply_group(r,net.ket_action)
    cases=[(;name="g$(g)q$(q)",right=g==0 ? r : phase_flip(other),sector=Z2Irrep(q),
        color=g==0 ? "blue" : "red",marker=q==0 ? "plus" : "circle") for g in 0:1 for q in 0:1]
    (;op=net.mpo,r,env,res,cases,actions=[net.ket_action,net.bra_action],basis="diagonal-Z2")
end
