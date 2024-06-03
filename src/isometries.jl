function CGmapping_from_AL(AL::TensorMap{T}, k0::Integer) where {T}
    D = dim(domain(AL))
    d = dims(codomain(AL))[2]

    D^2 < d^k0 || throw(ArgumentError("D^2 must be less than d^k0"))

    iDmat = diagm(ones(eltype(AL[]), D))

    indVecs = [[-(k0 + 1) -1 1], [[l -(l + 1) l + 1] for l in 1:k0-2]..., [(k0 - 1) -k0 -(k0 + 2)]]

    ψ_finite = ncon(fill(AL[], k0), indVecs, output=[-k0 - 2, -k0 - 1, -k0:1:-1...])

    ψ_finite = Matrix(reshape(ψ_finite, D^2, d^k0))
    F = qr(ψ_finite')
    V0k0 = Matrix(F.Q)'
    Pk0 = Matrix(F.R)'

    mpsInflatedR = zeros(eltype(AL[]), D^2, d, D^2)
    mpsInflatedL = zeros(eltype(AL[]), D^2, d, D^2)

    for l = 1:d
        mpsInflatedR[:, l, :] = kron(iDmat, AL[][:, l, :])
        mpsInflatedL[:, l, :] = kron(AL[][:, l, :], iDmat)
    end

    k = k0 + 1

    ψ_finite = ncon([mpsInflatedL, ψ_finite], [[-1, -3, 1], [1, -2]])
    _, r = qr(reshape(ψ_finite, D^2, d^k)')
    P = r'

    @tensor PmpsL[-1, -2, -3] := mpsInflatedL[-1, -3, 1] * Pk0[1, -2]
    PmpsL = reshape(PmpsL, D^2, (D^2) * d)
    @tensor PmpsR[-1, -2, -3] := mpsInflatedR[1, -2, -1] * Pk0[1, -3]
    PmpsR = reshape(PmpsR, D^2, (D^2) * d)

    L = P \ PmpsL
    R = P \ PmpsR

    return Matrix(V0k0), L, R
end

function good_ground_state(H::MPOHamiltonian{T},D::Integer) where {T}
    A = TensorMap(rand, ComplexF64, ℂ^D* ℂ^2, ℂ^D)
    state = InfiniteMPS([A])

    # groundstate, _, _= find_groundstate(
    #     state, H, VUMPS(;tol_galerkin=1e-4) 
    # )

    groundstate, _, _= find_groundstate(
        state, H, GradientGrassmann(;tol=1e-12) 
    )
    return groundstate
end


function approx_ground_state(H::MPOHamiltonian{T},ψ_good::MPSKit.AbstractMPS,d::Integer, D::Integer) where {T}
    data = TensorMap(rand,ComplexF64,ℂ^D*ℂ^d,ℂ^D)
    ψ_0 =  InfiniteMPS([data]);
    operator = make_time_mpo(H,0.1,WII())
    ψ_approx,_ = approximate(ψ_0,(operator,ψ_good),VUMPS(;tol_galerkin=1e-10))
    return ψ_approx
end