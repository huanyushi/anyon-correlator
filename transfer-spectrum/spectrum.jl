struct SectorLabel
    g::Int
    q::Int
    name::String
    function SectorLabel(g::Int,q::Int)
        g in (0,1) && q in (0,1) || throw(ArgumentError("g and q must be 0 or 1"))
        new(g,q,(("1","e"),("m","epsilon"))[g+1][q+1])
    end
end
function sector_backgrounds(pair,label)
    return (left=pair.r0,right=label.g==0 ? pair.r0 : pair.r1,
        leftenv=pair.env0,rightenv=label.g==0 ? pair.env0 : pair.env1,sector=Z2Irrep(label.q))
end

# Compatibility boundary: MPSKit 0.13.12, quasiparticleexcitation.jl lines
# 193-265 (statmech :LM and the public real-momentum Multiline wrapper),
# 282-294 and 374-389 (the SAME normalized effective action for residuals).
# Never call the InfiniteMPO/InfiniteQP excitations overload: it contains a
# self-reference at line 217 in this installed version.
struct SectorPoint
    g::Int
    q::Int
    k::Float64
    branch::Int
    ratio_real::Union{Nothing,Float64}
    ratio_imag::Union{Nothing,Float64}
    rho::Union{Nothing,Float64}
    phase::Union{Nothing,Float64}
    quasi_energy::Union{Nothing,Float64}
    projected_residual::Union{Nothing,Float64}
    variational_dim::Int
    status::String
    state
end
function mode_status(r,res,dim;zero_proven=false,budget=false,tol=1e-10,resolution=1e-12)
    budget && return "budget"
    dim==0 && return "empty_sector"
    (!isfinite(r) || !isfinite(res) || res>tol) && return "not_converged"
    zero_proven && return "zero_mode"
    abs(r)<=resolution && return "below_resolution"
    return "converged"
end
mode_phase(g,r)=g==0 ? angle(r) : nothing
function qp_coordinates(v)
    reduce(vcat,[vec(copy(b)) for x in v.Xs for (_,b) in blocks(x)])
end
function qp_from_coordinates(proto,coordinates)
    v=zero(proto);offset=0
    for x in v.Xs,(_,b) in blocks(x)
        n=length(b);b[:]=coordinates[offset+1:offset+n];offset+=n
    end
    offset==length(coordinates) || throw(DimensionMismatch("QP coordinates"))
    return v
end
function effective_action(op,l,r,proto;tol=1e-10)
    # l and r are the MPS backgrounds to the left/right of the defect B,
    # not the left/right eigenvectors of the non-Hermitian row transfer MPO.
    le=environments(l,op;tol); re=l===r ? le : environments(r,op;tol)
    E=MPSKit.effective_excitation_renormalization_energy(op,proto,le,re)
    eff=MPSKit.EffectiveExcitationHamiltonian(op,le,re,E)
    action=v->Base.invokelatest(eff,v;tol,maxiter=100,verbosity=0)
    return action,le,re
end
# Runtime boundaries keep Julia 1.12 from recursively inferring the whole
# eigensolver/environment graph into callers (a 600s top-level inference
# timeout was reproduced). Numerical algorithms and tolerance are unchanged.
function checked_spectrum_call(f)
    buffer=IOBuffer()
    result=with_logger(SimpleLogger(buffer,Logging.Warn)) do
        f()
    end
    messages=String(take!(buffer))
    failed=("GBL failed to converge","GBR failed to converge","excitation failed to converge",
        "linsolve stopped without converging","eigensolve stopped without converging","fixed point not converged after")
    any(s->occursin(s,messages),failed) && error("unconverged transfer excitation environment or eigensolve:\n"*messages)
    isempty(messages) || print(stderr,messages)
    return result
end
solve_modes(op,l,r,q,k;num=1,options=(;))=
    checked_spectrum_call(()->Base.invokelatest(_solve_modes,op,l,r,q,k;num,options))
function _solve_modes(op,l,r,q,k;num=1,options=(;))
    pkgversion(MPSKit)==v"0.13.12" || error("excitation adapter requires revalidation for this MPSKit version")
    tol=get(options,:tol,1e-10); resolution=get(options,:resolution,1e-12)
    reference=get(options,:reference,false)
    proto=Base.invokelatest(MPSKit.LeftGaugedQP,randn,l,r;sector=Z2Irrep(q),momentum=Float64(k))
    dimvar=sum(dim,proto.Xs)
    dimvar==0 && return (;ratios=ComplexF64[],residuals=Float64[],states=Any[],variational_dim=0,zero_proven=false)
    num>0 || throw(ArgumentError("num must be positive")); count=min(num,dimvar)
    action,le,re=Base.invokelatest(effective_action,op,l,r,proto;tol)
    # This is the full infinite-MPS single-QP parameter space, not a finite
    # circumference. Building every column retains degenerate multiplicities.
    limit=get(options,:dense_limit,2048)
    dimvar<=limit || throw(ArgumentError("QP dimension $dimvar exceeds explicit-matrix budget $limit"))
    projected=Matrix{ComplexF64}(undef,dimvar,dimvar)
    for j in 1:dimvar
        time()<get(options,:deadline,Inf) || error("spectrum matrix time budget reached")
        coordinates=zeros(ComplexF64,dimvar);coordinates[j]=1
        projected[:,j]=qp_coordinates(action(qp_from_coordinates(proto,coordinates)))
    end
    small_action=norm(projected)<=resolution
    # num controls how many eigenpairs are returned, not the matrix dimension.
    f=eigen(projected); order=sortperm(abs.(f.values);rev=true)[1:count]
    vals=ComplexF64.(f.values[order])
    states=[normalize!(qp_from_coordinates(proto,f.vectors[:,j])) for j in order]
    residuals=[Float64(norm(action(v)-value*v)) for (value,v) in zip(vals,states)]
    return (;ratios=vals,residuals,states,variational_dim=dimvar,zero_proven=reference && small_action)
end
function solve_sector(net,pair,label;momenta,num=1,options=(;))
    pair.diagnostics.converged || error("unconverged boundary pair")
    bg=sector_backgrounds(pair,label);points=SectorPoint[]
    tol=get(options,:tol,1e-10);resolution=get(options,:resolution,1e-12)
    residual_tol=get(options,:residual_tol,tol)
    deadline=get(options,:deadline,Inf)
    for k in momenta
        if time()>=deadline
            push!(points,SectorPoint(label.g,label.q,k,1,nothing,nothing,nothing,nothing,nothing,nothing,0,"budget",nothing));break
        end
        result=solve_modes(net.mpo,bg.left,bg.right,label.q,k;num,options)
        if result.variational_dim==0
            push!(points,SectorPoint(label.g,label.q,k,1,nothing,nothing,nothing,nothing,nothing,nothing,0,"empty_sector",nothing))
        end
        for (branch,(r,res,state)) in enumerate(zip(result.ratios,result.residuals,result.states))
            status=mode_status(r,res,result.variational_dim;zero_proven=result.zero_proven,tol=residual_tol,resolution)
            resolved=status=="converged" && abs(r)>resolution
            push!(points,SectorPoint(label.g,label.q,k,branch,real(r),imag(r),abs(r),resolved ? mode_phase(label.g,r) : nothing,
                resolved ? -log(abs(r)) : nothing,res,result.variational_dim,status,state))
        end
    end
    return points
end
