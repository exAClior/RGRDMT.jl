function load_Hamiltonian(filename::String)
    vars = matread(filename)
    return vars["H"], vars["Eexact"]
end

function load_MPS(filename::String, D::Integer)
    d = 2
    vars = matread(filename)

    upperBdFromMPS = vars["upperBdFromMPS"]
    MPSs = vars["MPS"] # why is there no imaginary part of MPS at all?

    MPSDn = MPSs[D]

    AL = TensorMap(Array(cat(MPSDn["AL"]..., dims=1)), ℂ^D * ℂ^d, ℂ^D)
    AR = TensorMap(Array(cat(MPSDn["AR"]..., dims=2)), ℂ^D * ℂ^d, ℂ^D)
    AC = TensorMap(Array(cat(MPSDn["AC"]..., dims=1)), ℂ^D * ℂ^d, ℂ^D)
    C = TensorMap(MPSDn["C"], ℂ^D, ℂ^D)
    ψ = InfiniteMPS([AL], [AR], [C], [AC])

    return ψ, upperBdFromMPS[D]
end