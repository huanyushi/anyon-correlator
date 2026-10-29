"Four-spin diagonal PEPS, with exp[(betax X + betaz Z)/4] on each edge."
function filtered_tensor(betax=.45,betaz=.55)
    # In the copy-V / XOR-P basis, a,b,c,d are the original vertex variables.
    # Split each copy V between the two selected checkerboard plaquettes.
    # Contracting four P and four half-vertex copy tensors gives this block K.
    # Equal zero-field bond weights contribute only an overall scalar.
    K=zeros(16,2,2,2,2)
    for a in 0:1,b in 0:1,c in 0:1,d in 0:1
        p=xor(a,b)+2xor(b,c)+4xor(c,d)+8xor(d,a)
        K[p+1,a+1,b+1,c+1,d+1]=1.
    end
    # Put the virtual X⊗4 symmetry into the Z-charge basis; HᵀH=I.
    A=K
    H=[1. 1.;1. -1.]/sqrt(2)
    for axis in 2:5;A=mode_array(A,H,axis);end
    F=exp([betaz betax;betax -betaz]/4)
    A=mode_array(A,kron(F,F,F,F),1)
    A/norm(A)
end
