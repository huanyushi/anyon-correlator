using Test
include("audit_state.jl")
const A=FrozenADAudit.A
@testset "frozen audit resumes only identical input and settings" begin
    mktempdir("/tmp/opencode") do dir
        b=A.ParameterBasis(2);x=A.initial_parameters(b,1)
        state=joinpath(dir,"state.json")
        A.save_state(state,x,b;metadata=(;hx=0.,hz=0.1,origin="symmetric_random",seed=1,accepted=false))
        out=joinpath(dir,"audit");mkpath(out)
        @test_throws ErrorException FrozenADAudit.main(["--state=$state","--out=$out"])
        settings=joinpath(out,"settings.json");write(settings,"{\"input_sha256\":\"different\"}")
        @test_throws ErrorException FrozenADAudit.main(["--state=$state","--out=$out","--resume-audit=true"])
        @test read(settings,String)=="{\"input_sha256\":\"different\"}"
        @test !isfile(joinpath(out,"accepted.json"))
    end
end
