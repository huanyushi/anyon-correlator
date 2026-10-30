"Orthonormal array bases for even C4 vertices and symmetric charge-preserving edges."
struct ParameterBasis
    D::Int
    charges::Vector{Int}
    V::Matrix{Float64}
    P::Matrix{Float64}
end
function ParameterBasis(D::Int;n_even=cld(D,2))
    D>=2 && 0<n_even<D || throw(ArgumentError("both virtual charges need positive multiplicity"))
    charges=vcat(zeros(Int,n_even),ones(Int,D-n_even))
    columns=Vector{Float64}[]; seen=Set{NTuple{4,Int}}()
    for i in CartesianIndices((D,D,D,D))
        t=Tuple(i)
        iseven(sum(charges[j] for j in t)) && !(t in seen) || continue
        orbit=unique([t,(t[2],t[3],t[4],t[1]),(t[3],t[4],t[1],t[2]),(t[4],t[1],t[2],t[3])])
        a=zeros(D,D,D,D)
        for j in orbit; a[j...]=inv(sqrt(length(orbit)));push!(seen,j);end
        push!(columns,vec(a))
    end
    vb=hcat(columns...);empty!(columns)
    for b in 1:D,a in 1:b,p in 1:2
        charges[a]==charges[b] || continue
        t=zeros(2,D,D);t[p,a,b]=t[p,b,a]=a==b ? 1. : inv(sqrt(2.))
        push!(columns,vec(t))
    end
    ParameterBasis(D,charges,vb,hcat(columns...))
end
parameter_count(b::ParameterBasis)=size(b.V,2)+size(b.P,2)+b.D
function parameter_pullback(x,b::ParameterBasis,g)
    y,back=Zygote.pullback(z->split_parameters(z,b),x)
    length(g)==length(y) || throw(DimensionMismatch("legacy gradient size"))
    result=Vector{Float64}(only(back(g)))
    all(isfinite,result) || throw(Legacy.GradientFailure("nonfinite shared gradient"))
    return result
end
function tensors(x,b::ParameterBasis)
    nv=size(b.V,2);np=size(b.P,2)
    Zygote.ignore() do
        length(x)==parameter_count(b) && all(isfinite,x) || throw(ArgumentError("invalid AD parameters"))
        norm(x[1:nv])>0 && norm(x[nv+1:nv+np])>0 || throw(ArgumentError("zero tensor norm"))
    end
    v=b.V*x[1:nv];p=b.P*x[nv+1:nv+np]
    u=x[nv+np+1:end];l=exp.(u.-sum(u)/b.D)
    Zygote.ignore() do
        all(isfinite,l) && all(>(0),l) || throw(ArgumentError("unsupported weight range"))
    end
    (;V=reshape(v/norm(v),ntuple(_->b.D,4)),P=reshape(p/norm(p),2,b.D,b.D),lambda=l)
end
function split_parameters(x,b::ParameterBasis)
    st=tensors(x,b)
    # Legacy AD splits each full connection weight into two square roots.
    # Repeat represents identical physical parameters, not independent sites.
    u=x[end-b.D+1:end]
    vcat(vec(st.V),vec(st.P),vec(st.P),repeat(u,4))
end
function initial_parameters(b::ParameterBasis,seed::Int;exact=false)
    exact || return vcat(randn(MersenneTwister(seed),size(b.V,2)+size(b.P,2)),zeros(b.D))
    # Embed the two-dimensional toric-code state into one copy of each charge.
    active=(1,findfirst(==(1),b.charges));v=zeros(ntuple(_->b.D,4));p=zeros(2,b.D,b.D)
    for t in Iterators.product(ntuple(_->active,4)...)
        iseven(sum(b.charges[j] for j in t)) && (v[t...]=1.)
    end
    for a in active,s in 0:1;p[s+1,a,a]=(-1.)^(s*b.charges[a]);end
    vcat(b.V'*vec(v),b.P'*vec(p),zeros(b.D))
end

# Portable parameters and tensors. Retain the existing format for saved states.
const STATE_FORMAT="one-site-z2-ad-v1"
function save_state(path,x,b::ParameterBasis;metadata)
    ispath(path) && throw(ArgumentError("refusing to overwrite $path"))
    st=tensors(x,b)
    data=Dict("format"=>STATE_FORMAT,"algorithm"=>"constrained_energy_ad",
        "D"=>b.D,"charges"=>b.charges,"parameters"=>collect(x),
        "V"=>vec(st.V),"P"=>vec(st.P),"lambda"=>st.lambda,
        "gauge"=>"balanced_square_root_connection_weights",
        "metadata"=>metadata)
    mkpath(dirname(abspath(path)))
    tmp,io=mktemp(dirname(abspath(path)))
    try
        JSON.print(io,data);flush(io);close(io)
        mv(tmp,path;force=false)
    finally
        isopen(io) && close(io)
        isfile(tmp) && rm(tmp)
    end
    return path
end
function load_state(path)
    data=JSON.parsefile(path)
    get(data,"format",nothing)==STATE_FORMAT &&
        get(data,"algorithm",nothing)=="constrained_energy_ad" ||
        throw(ArgumentError("not a supported constrained AD state"))
    b=ParameterBasis(Int(data["D"]);n_even=count(==(0),data["charges"]))
    data["charges"]==b.charges || throw(ArgumentError("charge ordering mismatch"))
    x=Float64.(data["parameters"]);st=tensors(x,b)
    data["gauge"]=="balanced_square_root_connection_weights" || throw(ArgumentError("unknown gauge"))
    for (key,expected) in (("V",vec(st.V)),("P",vec(st.P)),("lambda",st.lambda))
        actual=Float64.(data[key])
        length(actual)==length(expected) && isapprox(actual,expected;rtol=1e-13,atol=1e-14) ||
            throw(ArgumentError("stored $key disagrees with parameters"))
    end
    data["metadata"] isa AbstractDict || throw(ArgumentError("metadata must be an object"))
    return x,b,data["metadata"]
end
