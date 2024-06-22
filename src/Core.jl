function one_step_approx_dual(h::AbstractMatrix{V}, n::Integer, optimizer=SCS.Optimizer) where {V}

    if all(isreal, h)
        h = real(h)
    else
        throw(ArgumentError("Hamiltonian must be real"))
    end

    model = Model(optimizer)
    set_string_names_on_creation(model, false)
    d = 2 # spin physical dimension

    α = @variable(model, [1:d^(n-1), 1:d^(n-1)] in SymmetricMatrixSpace())

    @variable(model, ϵ)

    @expression(model, expr, kron(h, sparse(I, d^(n - 2), d^(n - 2))) + kron(sparse(I, d, d), α) - kron(α, sparse(I, d, d)) - ϵ .* sparse(I, d^n, d^n))

    @constraint(model, expr >= 0, PSDCone())
    @objective(model, Max, ϵ)

    optimize!(model)

    ElocTI = value(ϵ)
    expr = value.(expr)
    minEig, _, _ = eigsolve(expr, rand(eltype(expr), size(expr, 1)), 1, :SR)
    minEig = minEig[1]
    ElocTIRig = value(ϵ) + minEig
    return ElocTIRig, ElocTI
end

function one_step_approx(h::AbstractMatrix{V}, n::Integer, optimizer=SCS.Optimizer) where {V}
    model = Model(optimizer)
    set_string_names_on_creation(model, false)
    d = 2 # spin physical dimension
    ρs = [@variable(model, [1:d^sites, 1:d^sites] in HermitianPSDCone()) for sites in 2:n]

    @constraint(model, tr(ρs[1]) == 1.0)

    for ii in 3:n
        @constraint(model, ρs[ii-2] .== partialtrace(ρs[ii-1], 1, d * ones(Int64, ii)))
        @constraint(model, ρs[ii-2] .== partialtrace(ρs[ii-1], ii, d * ones(Int64, ii)))
    end

    @objective(model, Min, real(tr(h * ρs[1])))

    optimize!(model)
    return model
end

function two_step_approx(h::AbstractMatrix{V}, D::Integer, n::Integer,
    W2::AbstractMatrix{T}, L2::AbstractMatrix{T}, R2::AbstractMatrix{T},
    optimizer=SCS.Optimizer) where {V,T}

    model = Model(optimizer)
    set_string_names_on_creation(model, false)
    d = 2 # spin physical dimension
    k0 = floor(2 * log(D) / log(d)) + 1
    ρ3 = @variable(model, [1:d^(k0+1), 1:d^(k0+1)] in HermitianPSDCone())

    idmat = diagm(ones(eltype(h), d))

    ωs = [@variable(model, [1:d^2*D^2, 1:d^2*D^2] in HermitianPSDCone()) for _ in (k0+2):n]

    @constraint(model, tr(ρ3) == 1.0)
    @constraint(model,
        partialtrace(ρ3, 1, d * ones(Int64, k0 + 1)) .== partialtrace(ρ3, k0 + 1, d * ones(Int64, k0 + 1)))

    @constraint(model,
        kron(W2, idmat) * ρ3 * kron(W2, idmat)' .== partialtrace(ωs[1], 1, [d, D^2, d])
    )

    @constraint(
        model,
        kron(idmat, W2) * ρ3 * kron(idmat, W2)' .== partialtrace(ωs[1], 3, [d, D^2, d])
    )


    for ii in 2:n-(k0+1)
        @constraint(model, kron(L2, idmat) * ωs[ii-1] * kron(L2, idmat)' .== partialtrace(ωs[ii], 1, [d, D^2, d]))
        @constraint(model, kron(idmat, R2) * ωs[ii-1] * kron(idmat, R2)' .== partialtrace(ωs[ii], 3, [d, D^2, d]))
    end

    @objective(model, Min, real(tr(kron(h, diagm(ones(eltype(h), d^(k0 - 1)))) * ρ3)))

    optimize!(model)
    return model
end

