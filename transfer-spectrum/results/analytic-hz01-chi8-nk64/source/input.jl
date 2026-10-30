# Adapter from the research project's V/P/lambda representation to a blocked A.
transfer_arrays(st)=merge(dense_arrays(st),(;charges=[0,1]))
transfer_arrays(st::NamedTuple)=st
function build_transfer(st;direction=:down)
    direction==:down || throw(ArgumentError("only the confirmed downward transfer is defined"))
    a=transfer_arrays(st); D=length(a.lambda)
    ne=count(==(0),a.charges)
    a.charges==vcat(zeros(Int,ne),ones(Int,D-ne)) && 0<ne<D || throw(ArgumentError("invalid virtual charge order"))
    size(a.V)==(D,D,D,D) && size(a.P)==(2,D,D) || throw(DimensionMismatch("transfer input shape"))
    A=zeros(4,D,D,D,D)
    for p in 0:1,q in 0:1,n in 1:D,e in 1:D,s in 1:D,w in 1:D,b in 1:D,c in 1:D
        A[p+2q+1,n,e,s,w]+=a.V[c,b,s,w]*a.P[p+1,b,e]*a.P[q+1,c,n]*a.lambda[b]*a.lambda[e]*a.lambda[c]*a.lambda[n]
    end
    h=bytes2hex(sha256(reinterpret(UInt8,[vec(a.V);vec(a.P);a.lambda])))
    build_transfer_from_tensor(A,a.charges;input_hash=h)
end

# Zero-field label calibration, used by the research adapter and its tests.
function string_dictionary()
    p=dense_arrays(exact_reference_state()).P
    X=[0. 1;1 0];Z=[1. 0;0 -1];I2=Matrix{Float64}(I,2,2)
    virtual_error=norm(mode_array(p,X,1)-mode_array(p,Z,2))
    Ze=kron(I2,Z); Zw=kron(Z,I2); Xe=kron(I2,X); Xw=kron(X,I2)
    electric=Ze*Zw; magnetic=Xe*Xw
    e_end=max(norm(electric*Xe+Xe*electric),norm(electric*Xw+Xw*electric))
    e_bulk=norm(electric*(Xe*Xw)-(Xe*Xw)*electric)
    m_end=max(norm(magnetic*Ze+Ze*magnetic),norm(magnetic*Zw+Zw*magnetic))
    m_bulk=norm(magnetic*(Ze*Zw)-(Ze*Zw)*magnetic)
    validated=maximum((virtual_error,e_end,e_bulk,m_end,m_bulk))<=1e-12
    return (;validated,virtual_X_equals_Z=virtual_error,e_endpoint_anticommutator=e_end,
        e_interior_commutator=e_bulk,m_endpoint_anticommutator=m_end,m_interior_commutator=m_bulk)
end
