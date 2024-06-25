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

function two_step_approx(h::AbstractMatrix{V}, D::Integer, n::Integer,
    W2::AbstractMatrix{T}, L2::AbstractMatrix{T}, R2::AbstractMatrix{T},
    optimizer=SCS.Optimizer) where {V,T}

    model = Model(optimizer)
    set_string_names_on_creation(model, false)
    d = 2 # spin physical dimension
    k0 = Int(floor(2 * log(D) / log(d)) + 1)

    ρ = @variable(model, [1:d^(k0+1), 1:d^(k0+1)] in SymmetricMatrixSpace())

    ωs = [@variable(model, [1:d^2*D^2, 1:d^2*D^2] in SymmetricMatrixSpace()) for _ in k0:n-2]

    @constraint(model, ρ in PSDCone())
    for ii in 1:n-(k0+3)
        @constraint(model, ωs[ii] in PSDCone())
    end

    @constraint(model,
        partialtrace(ρ, 1, d * ones(Int64, k0 + 1)) .== partialtrace(ρ, k0 + 1, d * ones(Int64, k0 + 1)))

    @constraint(model, tr(ρ) == 1.0)

    sp_eyed = sparse(I, d, d)
    W2xI = kron(W2, sp_eyed)
    IxW2 = kron(sp_eyed, W2)

    @constraint(model,
        W2xI * ρ * W2xI' .== partialtrace(ωs[1], 1, [d, D^2, d])
    )

    @constraint(
        model,
        IxW2 * ρ * IxW2' .== partialtrace(ωs[1], 3, [d, D^2, d])
    )

    if k0 <= n -2
        LxI = kron(L2, sp_eyed) 
        IxR = kron(sp_eyed, R2)
        for ii in 2:n-(k0+1)
            @constraint(model, LxI * ωs[ii-1] * LxI' .== partialtrace(ωs[ii], 1, [d, D^2, d]))
            @constraint(model, IxR * ωs[ii-1] * IxR' .== partialtrace(ωs[ii], 3, [d, D^2, d]))
        end
    end


    @objective(model, Min, real(tr(kron(h, sparse(I, d^(k0 - 1), d^(k0 - 1))) * ρ)))

    optimize!(model)
    return objective_value(model)
end

