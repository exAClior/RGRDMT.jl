function main_dual(h::AbstractMatrix{V}, n_rng::AbstractRange,
    E_exact::Float64, efilename::String, nfilename::String,
    optimizer=SCS.Optimizer) where {V}

    vals = Float64[]
    for n in n_rng
        println("Working on n = $n")
        lower_bound, _ = one_step_approx_dual(h, n, optimizer)
        push!(vals, lower_bound)
        ΔELTI = real.(E_exact .- vals)
        writedlm(efilename, ΔELTI, ',')
        writedlm(nfilename, n_rng, ',')
    end
end

function main2(h::AbstractMatrix{T}, D::Integer,
    E_exact::Float64, n_rng::AbstractVector,
    W2::AbstractMatrix{T}, L2::AbstractMatrix{T}, R2::AbstractMatrix{T},
    efilename::String, nfilename::String,
    optimizer=MosekTools.Optimizer) where {T}

    vals = Float64[]
    for n in n_rng
        println("Working on n = $n")
        val = two_step_approx(h, D, n, W2, L2, R2, optimizer)
        push!(vals, val)
    end

    ΔErlxD = real.(E_exact .- vals)
    writedlm(efilename, ΔErlxD, ',')
    writedlm(nfilename, n_rng, ',')
end