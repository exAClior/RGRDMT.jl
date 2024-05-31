using Test, RGRDMT, Plots, Yao
using MPSKit, MPSKitModels
using TensorKit
using SCS, MosekTools, Dualization
using Random
using DelimitedFiles

function main(h::AbstractMatrix{V}, n_rng::UnitRange{Int}, E_exact::Float64, efilename::String, nfilename::String) where {V}
    vals = Float64[]
    for n in n_rng
        problem = one_step_approx(h, n, dual_optimizer(SCS.Optimizer))
        push!(vals, problem.optval)
    end
    ΔELTI = real.(E_exact .- vals)
    writedlm(efilename, ΔELTI, ',')
    writedlm(nfilename, n_rng, ',')
end


# Recreation of Fig2 (a)
J = 1.0
g = 1 / 2.0 # critical TFI
# E_exact = - 2/π *(1+g) * ellipe(4*g/(1+2*g+g^2))
H = transverse_field_ising(; J=J, g=g)
h = mat((kron(Z, Z)) + g / 2.0 * kron(X, I2) + g / 2.0 * kron(I2, X))
ψ, _, _ = mps_state(H, 2, 6)
E_exact = real(expectation_value(ψ, H))

main(h, 3:7, E_exact[1], "elti.csv", "n_exact.csv")

function main2(D::Integer, E_exact::Float64, n_rng::UnitRange{Int},A::Array{Complex{Float64},3})
    vals = Float64[]
    for n in n_rng
        println("Working on n = $n")
        problem = two_step_approx(h, D, n, A, dual_optimizer(Mosek.Optimizer)) # which optimizer did they use?
        push!(vals, problem.optval)
    end

    ΔErlxD = real.(E_exact .- vals)
    writedlm("erlxd$D.csv", ΔErlxD, ',')
    writedlm("nd$D.csv", n_rng, ',')
end

D = 3
d = 2
ψ, _, _ = mps_state(H, d, D)
A = ψ.AR[]
A = Array(reshape(A[],D,d,D)) 


Random.seed!(1234)
# A = rand_unitary(d * D)
# A = A[1:D, :]
# A * A'
# A = reshape(A, D, d, D)

main2(D, E_exact[1], 5:30,A)

function my_plot(eng_filenames::Vector{String}, n_filenames::Vector{String})
    plt = Plots.plot(
        xlabel="log2(n)",
        ylabel="log10(ΔE)",
        ylims=(1e-6, 1e-1),
        yticks=[1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1],
        yscale=:log10,
        xlims=(2, 50),
        # xticks = [2, 4,6,8,10,20,60,100,180],
        xscale=:log2,
        # xtickfont = font(20, "Courier")
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

# δ = 1.0
# H = heisenberg_XXZ(; Delta=δ , spin=1 // 2);
# h = mat((kron(X, X) + kron(Y, Y) + δ*kron(Z, Z)) / 4);
# ψ = mps_state(H, d, D)

# A = rand_unitary(d*D);
# A = A[1:D,:];
# A = reshape(A,D,d,D);
# typeof(A)

# n =  20
# res = main(h,D,n,A)

# @show res.optval