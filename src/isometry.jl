"""
    CGmapping_from_AL(AL::AbstractTensorMap, k0::Integer, n::Integer)

Constructs a CG mapping from an abstract tensor map `AL`.

# Arguments
- `AL::AbstractTensorMap`: The abstract tensor map.
- `k0::Integer`: The value of k0.
- `n::Integer`: The value of n.

# Returns
- If `n == k0`, returns a tuple `(V0, V1, V2)` where:
    - `V0`: A matrix of size D^2 x (d^k0) representing the finite part of the CG mapping.
    - `V1`: A sparse matrix of size D^2 x (D^2 * d) representing the left-infinite part of the CG mapping.
    - `V2`: A sparse matrix of size D^2 x (D^2 * d) representing the right-infinite part of the CG mapping.
- If `n != k0`, returns a tuple `(V0, L, R)` where:
    - `V0`: A matrix of size D^2 x (d^k0) representing the finite part of the CG mapping.
    - `L`: A matrix of size D^2 x (D^2 * d) representing the left-infinite part of the CG mapping.
    - `R`: A matrix of size D^2 x (D^2 * d) representing the right-infinite part of the CG mapping.
"""
function CGmapping_from_AL(AL::AbstractTensorMap, k0::Integer, n::Integer)
    D = dim(domain(AL))
    d = dims(codomain(AL))[2]

    D^2 <= d^k0 || throw(ArgumentError("D^2 must be less than d^k0"))

    sp_eyeD = sparse(I, D, D)

    mpsInflatedR = zeros(eltype(AL[]), D^2, d, D^2)
    mpsInflatedL = zeros(eltype(AL[]), D^2, d, D^2)

    for l = 1:d
        mpsInflatedR[:, l, :] = kron(sp_eyeD, AL[][:, l, :])
        mpsInflatedL[:, l, :] = kron(AL[][:, l, :], sp_eyeD)
    end

    indVecs = [[-(k0 + 1) -1 1], [[l -(l + 1) l + 1] for l in 1:k0-2]..., [(k0 - 1) -k0 -(k0 + 2)]]

    ψ_finite = ncon(fill(AL[], k0), indVecs, output=[-k0 - 2, -k0 - 1, -k0:1:-1...])

    V0 = Matrix(reshape(ψ_finite, D^2, d^k0))

    if n == k0
        return V0, sparse(one(eltype(V0)) * I, D^2, (D^2) * d), sparse(one(eltype(V0)) * I, D^2, D^2 * d)
    else
        L = reshape(permutedims(mpsInflatedL, (1, 3, 2)), D^2, (D^2) * d)
        R = reshape(permutedims(mpsInflatedR, (3, 2, 1)), D^2, (D^2) * d)
        return V0, L, R
    end
end

"""
    good_ground_state(H::MPOHamiltonian{T}, D::Integer) where {T}

Compute the ground state of a given MPOHamiltonian using the VUMPS algorithm.

# Arguments
- `H::MPOHamiltonian{T}`: The MPOHamiltonian representing the Hamiltonian of the system.
- `D::Integer`: The bond dimension of the MPS.

# Returns
- `groundstate`: The ground state of the system.
"""
function good_ground_state(H::MPOHamiltonian{T}, D::Integer) where {T}
    d = 2
    A = TensorMap(rand, ComplexF64, ℂ^D * ℂ^d, ℂ^D)
    A = normalize(A)
    state = InfiniteMPS([A])

    groundstate, _, _ = find_groundstate(
        state, H, VUMPS(; tol_galerkin=1e-14, maxiter=400, verbose=false)
    )

    return groundstate
end

"""
approx_ground_state(H::MPOHamiltonian{T}, ψ_good::InfiniteMPS, d::Integer, D::Integer) where {T}

Approximates the ground state of a given Hamiltonian using the variational uniform matrix product state (VUMPS) algorithm.

# Arguments
- `H::MPOHamiltonian{T}`: The Hamiltonian for which the ground state is to be approximated.
- `ψ_good::InfiniteMPS`: An initial guess for the ground state.
- `d::Integer`: The local Hilbert space dimension.
- `D::Integer`: The bond dimension.

# Returns
- `ψ_approx::InfiniteMPS`: The approximate ground state.

"""
function approx_ground_state(H::MPOHamiltonian{T}, ψ_good::InfiniteMPS, d::Integer, D::Integer) where {T}
    data = TensorMap(rand, ComplexF64, ℂ^D * ℂ^d, ℂ^D)
    ψ_0 = InfiniteMPS([data])
    operator = make_time_mpo(H, 0.1, WII())
    ψ_approx, _ = approximate(ψ_0, (operator, ψ_good), VUMPS(; tol_galerkin=1e-10, verbose=false))
    return ψ_approx
end