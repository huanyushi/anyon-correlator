include("core.jl")
options=Dict{String,String}()
for arg in ARGS
    startswith(arg,"--") && occursin('=',arg) || error("expected --key=value")
    key,value=split(arg[3:end],'=';limit=2);options[key]=value
end
OneSiteZ2AD.run(options)
