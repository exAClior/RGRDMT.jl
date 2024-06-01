using Test, RGRDMT, Plots, Yao
using MPSKit, MPSKitModels
using TensorKit
using SCS, MosekTools, Dualization
using Random
using DelimitedFiles

function main(h::AbstractMatrix{V}, n_rng::UnitRange{Int}, E_exact::Float64, efilename::String, nfilename::String) where {V}
    vals = Float64[]
    for n in n_rng
        println("Working on n = $n")
        problem = one_step_approx(h, n, dual_optimizer(SCS.Optimizer))
        push!(vals, problem.optval)
        ΔELTI = real.(E_exact .- vals)
        writedlm(efilename, ΔELTI, ',')
        writedlm(nfilename, n_rng, ',')
    end
end


# Recreation of Fig2 (a)
# J = 1.0
# g = 1 / 2.0 # critical TFI
# H = transverse_field_ising(; J=J, g=g)

# h = mat(-kron(Z, Z)/4.0 - g * (kron(X, I2)/2.0+ kron(I2, X)/2.0))

δ = 1.0
# H = heisenberg_XXZ(; Delta=δ , spin=1 // 2);
H = heisenberg_XXX(; spin=1 // 2);
h = mat((kron(X, X) + kron(Y, Y) + δ*kron(Z, Z)) / 4);

ψ, _, _ = mps_state(H, 2, 10);
E_exact = real(expectation_value(ψ, H))
E_exact = 0.25 - log(2)

main(h, 3:8, E_exact[1], "elti.csv", "n_exact.csv")

function main2(h::AbstractMatrix{T}, k0::Integer, D::Integer, E_exact::Float64, n_rng::UnitRange{Int}, W2::AbstractMatrix{T}, L2::AbstractMatrix{T}, R2::AbstractMatrix{T}) where {T}
    vals = Float64[]
    for n in n_rng
        println("Working on n = $n")
        problem = two_step_approx(h, k0, D, n, W2, L2, R2, Mosek.Optimizer) # which optimizer did they use?
        push!(vals, problem.optval)
    end

    ΔErlxD = real.(E_exact .- vals)
    writedlm("erlxd$D.csv", ΔErlxD, ',')
    writedlm("nd$D.csv", n_rng, ',')
end


D = 3
d = 2
ψ, _, _ = mps_state(H, d, D)

k0 = 4
A = ψ.AL[]
V0, L, R = CGmapping_from_AL(A, k0)


main2(h, k0, D, E_exact[1], 8:20, V0, L, R)


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


eng_filenames = ["elti.csv", "erlxd2.csv", "erlxd3.csv"]
n_filenames = ["n_exact.csv", "nd2.csv", "nd3.csv"]

cur_plt = my_plot(eng_filenames, n_filenames)

# ψ = mps_state(H, d, D)
