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
