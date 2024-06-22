using Test, RGRDMT, Yao
using MPSKit, MPSKitModels, TensorKit
using SCS, MosekTools, Dualization
using Random
using DelimitedFiles
using Plots
using JuMP

# @testset "Performance" begin
#     using BenchmarkTools
#     n = 6
#     optimizer = MosekTools.Optimizer
#     h = mat(-kron(Z, Z) / 4 - 1 / 4 * kron(X, I2))
#     problem = one_step_approx(h, n, dual_optimizer(optimizer)) # 7.74s 43.338MiB
#     problem = one_step_approx(h, n, optimizer) # 20s didn't terminate 66.720MiB
#     model = one_step_approx_jp(h, n, optimizer) # 13.53s 115212272 bytes
#     model = one_step_approx_jp(h, n, dual_optimizer(optimizer)) # 10.62s
# end

function main_jp(h::AbstractMatrix{V}, n_rng::UnitRange{Int},
    E_exact::Float64, efilename::String, nfilename::String,
    optimizer=SCS.Optimizer) where {V}

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

function main2_jp(h::AbstractMatrix{T}, k0::Integer, D::Integer,
    E_exact::Float64, n_rng::UnitRange{Int},
    W2::AbstractMatrix{T}, L2::AbstractMatrix{T}, R2::AbstractMatrix{T},
    optimizer=MosekTools.Optimizer) where {T}

    vals = Float64[]
    for n in n_rng
        println("Working on n = $n")
        model = two_step_approx(h, k0, D, n, W2, L2, R2, optimizer)
        push!(vals, objective_value(model))
    end

    ΔErlxD = real.(E_exact .- vals)
    writedlm("erlxd$D.csv", ΔErlxD, ',')
    writedlm("nd$D.csv", n_rng, ',')
end

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
    main_jp(h, 3:10, -1 / π, "data/etfi.csv", "data/ntfi.csv", dual_optimizer(MosekTools.Optimizer))

    eng_filenames = ["data/etfi.csv"]
    n_filenames = ["data/ntfi.csv"]
    # eng_filenames = ["elti.csv", "erlxd2.csv", "erlxd5.csv"]
    # n_filenames = ["n_exact.csv", "nd2.csv", "nd5.csv"]

    cur_plt = my_plot(eng_filenames, n_filenames)
end

booda()

function dooda()
    h = mat(-kron(X, X) - kron(Y, Y) + kron(Z, Z)) / 4
    main_jp(h, 3:10, 1 / 4 - log(2), "data/exxx.csv", "data/nxxx.csv", dual_optimizer(MosekTools.Optimizer))

    eng_filenames = ["data/exxx.csv"]
    n_filenames = ["data/nxxx.csv"]

    cur_plt = my_plot(eng_filenames, n_filenames)
end

dooda()