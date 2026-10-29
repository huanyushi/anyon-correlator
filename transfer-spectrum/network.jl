struct TransferNetwork
    A::Array{Float64,5}
    Q_raw::Array{Float64,4}
    Q::Array{Float64,4}
    mpo::InfiniteMPO
    fusion::TensorMap
    ket_action::TensorMap
    bra_action::TensorMap
    step::Int
    normalization::Float64
    input_hash::String
end
function mode_array(a,m,axis)
    order=[axis;[i for i in 1:ndims(a) if i!=axis]]; dims=size(a)[order]
    permutedims(reshape(m*reshape(permutedims(a,order),dims[1],:),size(m,1),dims[2:end]...),invperm(order))
end
"Ordinary down-row MPO from an already blocked physical,n,e,s,w tensor."
function build_transfer_from_tensor(A::Array{Float64,5},charges;input_hash=bytes2hex(sha256(reinterpret(UInt8,vec(A)))))
    D=length(charges);d=D^2;ne=count(==(0),charges)
    charges==vcat(zeros(Int,ne),ones(Int,D-ne)) && 0<ne<D || throw(ArgumentError("invalid virtual charge order"))
    size(A)[2:5]==(D,D,D,D) || throw(DimensionMismatch("blocked tensor shape"))
    @tensor raw[n,N,e,E,s,S,w,W] := A[p,n,e,s,w]*conj(A[p,N,E,S,W])
    Q_raw=reshape(raw,d,d,d,d)
    v=Z2Space(0=>ne,1=>D-ne); W=fuse(v⊗conj(v)); F=unitary(Float64,W,v⊗conj(v))
    f=reshape(convert(Array,F),d,d); Q=copy(Q_raw)
    for i in 1:4; Q=mode_array(Q,f,i); end
    signs=Matrix(Diagonal(Float64.(1 .-2charges)))
    z=TensorMap(signs,v←v)
    ket=F*(z⊗id(Float64,conj(v)))*F'
    zb=TensorMap(signs,conj(v)←conj(v))
    bra=F*(id(Float64,v)⊗zb)*F'
    scale=norm(Q); scale>0 || error("zero double layer")
    # Ordinary MPO: left,down <- up,right. No seam or twist is inserted.
    local_mpo=TensorMap(permutedims(Q/scale,(4,3,1,2)),W⊗W←W⊗W;tol=1e-12)
    TransferNetwork(A,Q_raw,Q,InfiniteMPO([local_mpo]),F,ket,bra,1,scale,input_hash)
end
"Tiny periodic-row algebra check, never the infinite-plane production spectrum."
function closed_row(net,width)
    size(net.Q,1)==4 || throw(ArgumentError("closed-row algebra checks require D=2"))
    width in (2,3) || throw(ArgumentError("algebra checks only support widths 2,3"))
    q=net.Q/net.normalization; d=4^width; t=zeros(d,d)
    tuples=collect(Iterators.product(ntuple(_->1:4,width)...))[:]
    for (i,s) in enumerate(tuples),(j,n) in enumerate(tuples),bonds in tuples
        value=1.0
        for x in 1:width; value*=q[n[x],bonds[x],s[x],bonds[mod1(x-1,width)]]; end
        t[i,j]+=value
    end
    return t
end
