# Read-only independent audit of the ALREADY measured 700-step state.
# Does not repeat optimization or CTMRG; original budget summary is unchanged.
const SU=normpath(joinpath(@__DIR__,"..",".."))
include(joinpath(SU,"core.jl"))
using .OneSiteZ2SU, SHA
const S=OneSiteZ2SU
const R=joinpath(SU,"results","zero-z2-seed5-final-audit")
st,md=load_state(joinpath(R,"latest.json"))
md["step"]==700 || error("wrong checkpoint")
events=[S.JSON.parse(s) for s in eachline(joinpath(R,"trajectory.jsonl"))]
measurements=[e["m"] for e in events if e["event"]=="measurement"]
sort([m["chi"] for m in measurements])==[8,16] || error("missing stored CTMRG measurements")
ms=[(;energy=Float64(m["energy"]),vals=Float64.(m["vals"]),residual=Float64(m["residual"]),
    converged=Bool(m["converged"]),chi=Int(m["chi"])) for m in measurements]
result=S.audit(ms,st,"finished")
original=S.JSON.parsefile(joinpath(R,"summary.json"))
out=joinpath(@__DIR__,"physics-audit-700.json")
ispath(out) && error("refusing to overwrite physics-audit-700.json")
S.write_json(out,(;record_type="independent read-only physics audit of saved state and measurements",
    audit_completed=true,accepted=result.passed,original_runner_status=original["status"],
    state_sha256=bytes2hex(sha256(read(joinpath(R,"latest.json")))),step=700,
    dt=md["dt"],seed=md["seed"],als_sweeps=md["als_sweeps"],measurements=ms,result,
    note="No new CTMRG, no optimization; local tensor residuals recomputed. Original metadata remains accepted=false."))
S.JSON.print(stdout,result);println();flush(stdout)
