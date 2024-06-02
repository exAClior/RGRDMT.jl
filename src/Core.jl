
function one_step_approx(h::AbstractMatrix{V}, n::Integer, optimizer=SCS.Optimizer) where {V}
    d = 2 # spin physical dimension
    ρs = [HermitianSemidefinite(d^ii, d^ii) for ii in 2:n]

    constraints = Constraint[
        tr(ρs[1])==1.0,
        [ρs[ii-2] == partialtrace(ρs[ii-1], 1, d * ones(Int64, ii)) for ii in 3:n]...,
        [ρs[ii-2] == partialtrace(ρs[ii-1], ii, d * ones(Int64, ii)) for ii in 3:n]...,
    ]
    problem = minimize(real(tr(h * ρs[1])), constraints)

    solve!(problem, optimizer; silent=true)
    return problem
end


function two_step_approx(h::AbstractMatrix{V}, k0::Integer, D::Integer, n::Integer, W2::AbstractMatrix{T}, L2::AbstractMatrix{T}, R2::AbstractMatrix{T}, optimizer=SCS.Optimizer) where {V,T}
    d = 2 # spin physical dimension

    D^2 < d^(k0) || throw(ArgumentError("D^2 must be greater than d^k0"))

    ρ3 = ComplexVariable(d^(k0 + 1), d^(k0 + 1))

    idmat = diagm(ones(eltype(h), d))

    ωs = [ComplexVariable(d^2 * D^2, d^2 * D^2) for _ in (k0+2):n]

    constraints = Constraint[
        tr(ρ3)==1.0,
        ρ3 ⪰ 0,
        [ω ⪰ 0 for ω in ωs]...,
        partialtrace(ρ3, 1, d * ones(Int64, k0 + 1))==partialtrace(ρ3, k0 + 1, d * ones(Int64, k0 + 1)),
        kron(W2, idmat)*ρ3*kron(W2, idmat)'==partialtrace(ωs[1], 1, [d, D^2, d]),
        kron(idmat, W2)*ρ3*kron(idmat, W2)'==partialtrace(ωs[1], 3, [d, D^2, d])
    ]

    for ii in 2:n-(k0+1)
        push!(constraints, kron(L2, idmat) * ωs[ii-1] * kron(L2, idmat)' == partialtrace(ωs[ii], 1, [d, D^2, d]))
        push!(constraints, kron(idmat, R2) * ωs[ii-1] * kron(idmat, R2)' == partialtrace(ωs[ii], 3, [d, D^2, d]))
    end

    problem = minimize(real(tr(kron(h, diagm(ones(eltype(h), d^(k0 - 1)))) * ρ3)), constraints)

    solve!(problem, optimizer; silent=true)
    return problem
end

