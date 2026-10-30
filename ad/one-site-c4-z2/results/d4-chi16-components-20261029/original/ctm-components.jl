include("/home/huanyushi/code/anyon-correlator/ad/one-site-z2/core.jl")
using .OneSiteZ2AD, LinearAlgebra, TensorKit, PEPSKit, JLD2
const A=OneSiteZ2AD
root="/home/huanyushi/code/anyon-correlator/"
x,b,md=A.load_state(root*"ad/one-site-z2/results/d4-exact-noise-cont1-20261029/final.json")
psi=A.Legacy.make_peps(A.split_parameters(x,b),b.D)
env=JLD2.load(root*"ad/one-site-z2/results/d4-chi16-trace-20261029/env-4000.jld2","env")
net=PEPSKit.InfiniteSquareNetwork(psi)
vals(t)=reduce(vcat,[collect(v) for (_,v) in blocks(PEPSKit.svd_vals(t))])
method=Symbol(get(ENV,"CTM_SVD","DefaultAlgorithm"))
alg=SimultaneousCTMRG(;tol=1e-10,maxiter=1,verbosity=0,decomposition_alg=(;alg=method))
println("SVD method=$method")
for step in 1:5
    newer,info=PEPSKit.ctmrg_iteration(net,env,alg)
    for name in (:corners,:edges)
        old=getproperty(env,name);new=getproperty(newer,name)
        errors=[norm(vals(a)-vals(b)) for (a,b) in zip(old,new)]
        idx=argmax(errors)
        println("step=$step component=$name index=$idx delta=$(errors[idx])")
        println("singular_values=",vals(new[idx]))
        println("changes=",vals(new[idx])-vals(old[idx]))
    end
    println("projector_singular_values=",[diag(convert(Array,s)) for s in info.S])
    println("metrics=",info.contraction_metrics);flush(stdout)
    global env=newer
end
