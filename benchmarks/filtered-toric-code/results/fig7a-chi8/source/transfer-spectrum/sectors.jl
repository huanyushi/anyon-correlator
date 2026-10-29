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
