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
# E_exact = 2/π
# h = mat(-kron(Z, Z)/4.0 - g * (kron(X, I2)/2.0+ kron(I2, X)/2.0))


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


function actual()

    h = mat((kron(X, X) + kron(Y, Y) + kron(Z, Z)) / 4)
    H = heisenberg_XXX(; spin=1 // 2)
    E_exact = 0.25 - log(2)
    main(h, 3:12, E_exact, "elti.csv", "n_exact.csv")


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

actual()

function my_plot(eng_filenames::Vector{String}, n_filenames::Vector{String})
    plt = Plots.plot(
        xlabel="log2(n)",
        ylabel="log10(ΔE)",
        ylims=(1e-6, 1e-1),
        # ylims=(1e-3, 1e-1),
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

eng_filenames = ["elti.csv"]
n_filenames = ["n_exact.csv"]
# eng_filenames = ["elti.csv", "erlxd2.csv", "erlxd5.csv"]
# n_filenames = ["n_exact.csv", "nd2.csv", "nd5.csv"]

cur_plt = my_plot(eng_filenames, n_filenames)


D = 120
A = TensorMap(rand, ComplexF64, ℂ^D * ℂ^2, ℂ^D)
state = InfiniteMPS([A])

XXXH = heisenberg_XXX(; spin=1 // 2)

optimize_steps = Float64[]
function finalize(iter, ψ, H, envs)
    push!(optimize_steps, real(expectation_value(ψ,H)[1]))
    return ψ, envs 
end
groundstate, cache,delta = find_groundstate(
    state, XXXH, VUMPS(;tol_galerkin=1e-5, verbose=false,finalize=finalize) 
)

E_exact = 1/4 - log(2)

δEs = optimize_steps .- E_exact

using Plots
Plots.plot(1:length(δEs), δEs, yscale=:log10, ylabel="log10(ΔE)", xlabel="iteration", label="VUMPS")
