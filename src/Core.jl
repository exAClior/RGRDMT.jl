function mps_state(H::MPOHamiltonian{T}, d::Integer, D::Integer) where {T}
    random_data = TensorMap(rand, ComplexF64, ℂ^D * ℂ^d, ℂ^D)
    state = InfiniteMPS([random_data])
    return find_groundstate(state, H, GradientGrassmann(; maxiter=50))
end


function one_step_approx_nativedual(h::AbstractMatrix{V}, n::Integer, optimizer=SCS.Optimizer) where {V}
    d = Int(sqrt(size(h,1)))   # spin physical dimension

    constraints = Constraint[
        tr(ρs[1]) == 1.0,
        [ρ ⪰ 0 for ρ in ρs]...,
        [ρs[ii-2] == partialtrace(ρs[ii-1], 1, d * ones(Int64, ii)) for ii in 3:n]...,
        [ρs[ii-2] == partialtrace(ρs[ii-1], ii, d * ones(Int64, ii)) for ii in 3:n]...,
    ]

    problem = minimize(real(tr(h * ρs[1])), constraints)

    solve!(problem, optimizer;silent=true)
    return problem
end


function one_step_approx(h::AbstractMatrix{V}, n::Integer, optimizer=SCS.Optimizer) where {V}
    d = 2 # spin physical dimension
    ρs = [HermitianSemidefinite(d^ii, d^ii) for ii in 2:n]

    constraints = Constraint[
        tr(ρs[1]) == 1.0,
        [ρs[ii-2] == partialtrace(ρs[ii-1], 1, d * ones(Int64, ii)) for ii in 3:n]...,
        [ρs[ii-2] == partialtrace(ρs[ii-1], ii, d * ones(Int64, ii)) for ii in 3:n]...,
    ]

    problem = minimize(real(tr(h * ρs[1])), constraints)

    solve!(problem, optimizer;silent=true)
    return problem
end

function two_step_approx(h::AbstractMatrix{V}, D::Integer, n::Integer, A::AbstractArray{V,3}, optimizer=SCS.Optimizer) where {V}
    d = 2 # spin physical dimension
    # 1d translational invariant hamiltonian local term

    ρ3 = ComplexVariable(d^3, d^3)

    @tensor W2[j, l, i, k] := A[j, i, a] * A[a, k, l]

    W2 = reshape(W2, D^2, d^2)

    R2 = reshape(A, D, D * d)
    L2 = transpose(reshape(A, D * d, D))

    iremain = Diagonal(ones(ComplexF64, d * D))
    i2mat = mat(I2)

    ωs = [ComplexVariable(d^2 * D^2, d^2 * D^2) for _ in 1:n-3]

    constraints = Constraint[
        tr(ρ3) == 1.0,
        ρ3 ⪰ 0,
        [ω ⪰ 0 for ω in ωs]...,
        partialtrace(ρ3, 1, d * ones(Int64, 3)) == partialtrace(ρ3, 3, d * ones(Int64, 3)),
        kron(W2, i2mat) * ρ3 * kron(W2, i2mat)' == partialtrace(ωs[1], 1, [d, D^2, d]),
        kron(i2mat, W2) * ρ3 * kron(i2mat, W2)' == partialtrace(ωs[1], 3, [d, D^2, d])
    ]

    for ii in 2:n-3
        push!(constraints, kron(L2, iremain) * ωs[ii-1] * kron(L2, iremain)' == partialtrace(ωs[ii], 1, [d, D^2, d]))
        push!(constraints, kron(iremain, R2) * ωs[ii-1] * kron(iremain, R2)' == partialtrace(ωs[ii], 3, [d, D^2, d]))
    end

    problem = minimize(real(tr(kron(h, mat(I2)) * ρ3)), constraints)

    solve!(problem, optimizer)
    return problem
end

