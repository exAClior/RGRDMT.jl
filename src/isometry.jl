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
        return V0, typeof(mpsInflatedL)[], typeof(mpsInflatedR)[]
    else
        L = reshape(permutedims(mpsInflatedL, (1, 3, 2)), D^2, (D^2) * d)
        R = reshape(permutedims(mpsInflatedR, (3, 2, 1)), D^2, (D^2) * d)
        return V0, L, R
    end
end

function good_ground_state(H::MPOHamiltonian{T}, D::Integer) where {T}
    d = 2
    A = TensorMap(rand, ComplexF64, ℂ^D * ℂ^d, ℂ^D)
    A = normalize(A)
    state = InfiniteMPS([A])

    groundstate, _, _ = find_groundstate(
        state, H, VUMPS(; tol_galerkin=1e-10, maxiter=200)
    )

    return groundstate
end

function approx_ground_state(H::MPOHamiltonian{T}, ψ_good::InfiniteMPS, d::Integer, D::Integer) where {T}
    data = TensorMap(rand, ComplexF64, ℂ^D * ℂ^d, ℂ^D)
    ψ_0 = InfiniteMPS([data])
    operator = make_time_mpo(H, 0.1, WII())
    ψ_approx, _ = approximate(ψ_0, (operator, ψ_good), VUMPS(; tol_galerkin=1e-10))
    return ψ_approx
end