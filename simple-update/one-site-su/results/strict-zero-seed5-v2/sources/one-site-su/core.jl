module OneSiteSU
using LinearAlgebra

# This is the actual successful zero-field implementation, not a reimplementation.
include("../split/core.jl")

"A spatially symmetric random start, not the analytic toric-code tensor."
function random_state(seed)
    st=SplitSU.initial_state(:random,seed;L=2)
    v=st.tensors[(:V,1,1)]
    v=(v+permutedims(v,(2,3,4,1))+permutedims(v,(3,4,1,2))+
       permutedims(v,(4,1,2,3)))/4
    p=st.tensors[(:E,1,1)]
    p=(p+permutedims(p,(1,3,2)))/2
    normalize!(v);normalize!(p)
    for n in keys(st.tensors)
        st.tensors[n]=copy(n[1]==:V ? v : p)
    end
    return st
end

"Direct call: original star HOSVD, plaquette ALS, and lambda restoration."
sweep!(st,dt,step;als_sweeps=12)=
    SplitSU.sweep!(st,dt,step;plaquette=:loop,als_sweeps)

"Raw representation differences, not a test for physical symmetry breaking."
function symmetry_defects(st)
    vs=[st.tensors[(:V,r,c)] for r in 1:st.L,c in 1:st.L]
    ps=[a for (n,a) in st.tensors if n[1]!=:V]
    ws=[st.weights[k] for k in sort(collect(keys(st.weights)))]
    relative(a,b)=norm(a-b)/max(norm(b),eps())
    v_shared=maximum(relative(a,first(vs)) for a in vs)
    p_shared=maximum(relative(a,st.tensors[(:E,1,1)]) for a in ps)
    lambda_shared=maximum(relative(w,first(ws)) for w in ws)
    v_rotation=maximum(relative(permutedims(a,(2,3,4,1)),a) for a in vs)
    p_reversal=maximum(relative(permutedims(a,(1,3,2)),a) for a in ps)
    return (;v_shared,p_shared,lambda_shared,v_rotation,p_reversal,
             maximum=max(v_shared,p_shared,lambda_shared,v_rotation,p_reversal))
end
end
