# Reproducible SOFTWARE calibration only; not a numerical PEPS ground-state run.
const TD=normpath(joinpath(@__DIR__,"..",".."))
include(joinpath(TD,"core.jl"))
using .ToricTransferSpectrum, TensorKit, LinearAlgebra, SHA
const T=ToricTransferSpectrum
const S=T.OneSiteZ2SU
function fixture(scale)
    W=Z2Space(0=>2,1=>2); b=Z2Space(0=>1)
    op=T.MPSKit.InfiniteMPO([TensorMap(reshape(scale*Matrix(Diagonal([5.,3,2,1])),1,4,4,1),b⊗W←W⊗b)])
    psi=T.MPSKit.InfiniteMPS([TensorMap(reshape(ComplexF64[1,0,0,0],1,4,1),b⊗W←b)])
    op,psi
end
out=joinpath(@__DIR__,"calibration.json")
ispath(out) && error("refusing to overwrite calibration.json")
records=Any[]
for scale in (1.,7.), k in (0.,pi/2), (q,num,expected) in ((0,1,[.6]),(1,2,[.4,.2]),(1,1,[.4]))
    op,psi=fixture(scale)
    result=T.solve_modes(op,psi,psi,q,k;num)
    err=maximum(abs.(result.ratios.-expected))
    row=(;scale,k,q,num,ratios_real=real.(result.ratios),ratios_imag=imag.(result.ratios),
        projected_residuals=result.residuals,max_absolute_error=err,variational_dim=result.variational_dim)
    push!(records,row);S.JSON.print(stdout,row);println();flush(stdout)
end
net=build_transfer(S.exact_reference_state());pair=solve_boundaries(net)
reference=Any[]
for g in 0:1,q in 0:1
    pts=solve_sector(net,pair,SectorLabel(g,q);momenta=[0.,pi/2],num=1,options=(;reference=true))
    append!(reference,T.point_data.(pts))
end
S.write_json(out,(;kind="software-control-and-exact-reference-only",compile_enabled=Base.JLOptions().compile_enabled,
    julia=string(VERSION),mpskit=string(pkgversion(T.MPSKit)),tensorkit=string(pkgversion(TensorKit)),
    tol=1e-10,resolution=1e-12,source_sha256=bytes2hex(sha256(read(joinpath(TD,"spectrum.jl")))),
    nonzero_controls=records,exact_boundary_diagnostics=pair.diagnostics,reference_points=reference))
println("CALIBRATION saved $out");flush(stdout)
