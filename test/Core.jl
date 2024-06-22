using Test, RGRDMT, Yao
using MPSKit, MPSKitModels, TensorKit
using SCS, MosekTools, Dualization
using Random
using DelimitedFiles
using Plots
using JuMP

@testset "Performance" begin
    using BenchmarkTools
    n = 6
    optimizer = MosekTools.Optimizer
    h = mat(-kron(Z, Z) / 4 - 1 / 4 * kron(X, I2))
    problem = one_step_approx(h, n, dual_optimizer(optimizer)) # 7.74s 43.338MiB
    problem = one_step_approx(h, n, optimizer) # 20s didn't terminate 66.720MiB
    model = one_step_approx_jp(h, n, optimizer) # 13.53s 115212272 bytes
    model = one_step_approx_jp(h, n, dual_optimizer(optimizer)) # 10.62s
end

function main(h::AbstractMatrix{V}, n_rng::UnitRange{Int}, E_exact::Float64, efilename::String, nfilename::String, optimizer=SCS.Optimizer) where {V}
    vals = Float64[]
    for n in n_rng
        println("Working on n = $n")
        problem = one_step_approx(h, n, dual_optimizer(optimizer))
        push!(vals, problem.optval)
        ΔELTI = real.(E_exact .- vals)
        writedlm(efilename, ΔELTI, ',')
        writedlm(nfilename, n_rng, ',')
    end
end


function main_jp(h::AbstractMatrix{V}, n_rng::UnitRange{Int}, E_exact::Float64, efilename::String, nfilename::String, optimizer=SCS.Optimizer) where {V}
    vals = Float64[]
    for n in n_rng
        println("Working on n = $n")
        model = one_step_approx_jp(h, n, optimizer)
        push!(vals, objective_value(model))
        ΔELTI = real.(E_exact .- vals)
        writedlm(efilename, ΔELTI, ',')
        writedlm(nfilename, n_rng, ',')
    end
end



function main2(h::AbstractMatrix{T}, k0::Integer, D::Integer, E_exact::Float64, n_rng::UnitRange{Int}, W2::AbstractMatrix{T}, L2::AbstractMatrix{T}, R2::AbstractMatrix{T}, optimizer=MosekTools.Optimizer) where {T}
    vals = Float64[]
    for n in n_rng
        println("Working on n = $n")
        problem = two_step_approx(h, k0, D, n, W2, L2, R2, optimizer)
        push!(vals, problem.optval)
    end

    ΔErlxD = real.(E_exact .- vals)
    writedlm("erlxd$D.csv", ΔErlxD, ',')
    writedlm("nd$D.csv", n_rng, ',')
end


function actual()

    h = mat((-kron(X, X) - kron(Y, Y) + kron(Z, Z)) / 4)
    H = heisenberg_XXX(; spin=1 // 2)
    E_exact = 0.25 - log(2)
    main(h, 3:12, E_exact, "exxx.csv", "nxxx.csv")


    # ψ = good_ground_state(H,100) 
    ψ = good_ground_state(H, 30)

    d = 2
    D = 2
    ψ_approx = approx_ground_state(H, ψ, d, D)
    # E_exact = real(expectation_value(ψ_approx, H))
    k0 = 3
    A = ψ_approx.AL[]
    V0, L, R = CGmapping_from_AL(A, k0)
    main2(h, k0, D, E_exact, 8:20, V0, L, R)

    D = 5
    ψ_approx = approx_ground_state(H, ψ, d, D)
    # E_exact = real(expectation_value(ψ_approx, H))
    k0 = 5
    A = ψ_approx.AL[]
    V0, L, R = CGmapping_from_AL(A, k0)
    main2(h, k0, D, E_exact, 8:30, V0, L, R)

    D = 4
    A = rand_unitary(d * D)
    A = A[1:D, :]
    A = reshape(A, D, d, D)
    A = TensorMap(A, ℂ^D * ℂ^d, ℂ^D)
    k0 = 5
    V0, L, R = CGmapping_from_AL(A, k0)
    main2(h, k0, D, E_exact, 8:30, V0, L, R)
end

# actual()

function my_plot(eng_filenames::Vector{String}, n_filenames::Vector{String})
    plt = Plots.plot(
        xlabel="log2(n)",
        ylabel="log10(ΔE)",
        ylims=(1e-6, 1e-1),
        yticks=[1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1],
        yscale=:log10,
        xlims=(2, 64),
        xscale=:log2,
    )
    for (efname, nfname) in zip(eng_filenames, n_filenames)
        ΔE = readdlm(efname, ',')
        n = readdlm(nfname, ',')
        Plots.plot!(plt, n, ΔE, label=efname[1:end-7], legend=:topright)
    end
    return plt
end


function booda()
    h = mat(-kron(Z, Z) / 4 - 1 / 4 * kron(X, I2))
    main_jp(h, 3:10, -1 / π, "data/etfi.csv", "data/ntfi.csv", MosekTools.Optimizer)
    main(h, 3:10, -1 / π, "data/etfi.csv", "data/ntfi.csv", MosekTools.Optimizer)

    eng_filenames = ["data/etfi.csv"]
    n_filenames = ["data/ntfi.csv"]
    # eng_filenames = ["elti.csv", "erlxd2.csv", "erlxd5.csv"]
    # n_filenames = ["n_exact.csv", "nd2.csv", "nd5.csv"]

    cur_plt = my_plot(eng_filenames, n_filenames)
end

booda()

function dooda()
    h = mat(-kron(X, X) - kron(Y, Y) + kron(Z, Z)) / 4
    main(h, 3:10, 1 / 4 - log(2), "data/exxx.csv", "data/nxxx.csv", MosekTools.Optimizer)

    eng_filenames = ["data/exxx.csv"]
    n_filenames = ["data/nxxx.csv"]

    cur_plt = my_plot(eng_filenames, n_filenames)
end

dooda()