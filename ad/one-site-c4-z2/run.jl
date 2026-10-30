include("core.jl")
options=Dict{String,String}()
for arg in ARGS
    startswith(arg,"--") && occursin('=',arg) || error("expected --key=value")
    key,value=split(arg[3:end],'=';limit=2);options[key]=value
end
get(options,"init","exact")=="exact" || error("random initial states use random.jl")
haskey(options,"seeds") && error("seed lists are only used by random.jl")
OneSiteZ2AD.run(options)
