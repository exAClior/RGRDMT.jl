"""
    load_Hamiltonian(filename::String)

Load the Hamiltonian matrix and the exact eigenvalues from a file.

This function is only used when loading ground state MPS from MATLAB computation using VUMPS implemented [here](https://github.com/Darkdragon84/IMPS_ML_tools)

# Arguments
- `filename::String`: The name of the file containing the Hamiltonian matrix.

# Returns
- `H`: The Hamiltonian matrix.
- `Eexact`: The exact eigenvalues.

"""
function load_Hamiltonian(filename::String)
    vars = matread(filename)
    return vars["H"], vars["Eexact"]
end

"""
    load_MPS(filename::String, D::Integer)

Load an MPS (Matrix Product State) from a file.

# Arguments
- `filename::String`: The path to the file containing the MPS data.
- `D::Integer`: The virtual bond dimension of the MPS.

# Returns
- `ψ`: The loaded MPS as an `InfiniteMPS` object.
- `upperBdFromMPS`: The upper bound of ground state energy achieved by the MPS.

"""
function load_MPS(filename::String, D::Integer)
    d = 2
    vars = matread(filename)

    upperBdFromMPS = vars["upperBdFromMPS"]
    MPSs = vars["MPS"] 

    MPSDn = MPSs[D]

    AL = TensorMap(Array(cat(MPSDn["AL"]..., dims=1)), ℂ^D * ℂ^d, ℂ^D)
    AR = TensorMap(Array(cat(MPSDn["AR"]..., dims=2)), ℂ^D * ℂ^d, ℂ^D)
    AC = TensorMap(Array(cat(MPSDn["AC"]..., dims=1)), ℂ^D * ℂ^d, ℂ^D)
    C = TensorMap(MPSDn["C"], ℂ^D, ℂ^D)
    ψ = InfiniteMPS([AL], [AR], [C], [AC])

    return ψ, upperBdFromMPS[D]
end