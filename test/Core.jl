using Test, RGRDMT, Plots, Yao
using MPSKit, MPSKitModels
using SCS, MosekTools
using Random

Random.seed!(1234)
d = 2

# Recreation of Fig2 (a)
J = 1.0
g = 1 / 2.0 # critical TFI
# E_exact = -g * (1.0 + J^2/4/g^2)
H = transverse_field_ising(; J=J, g=g)
h = mat((kron(Z, Z)) + g/2.0 * kron(X,I2) + g/2.0 * kron(I2,X));
ψ, _ , _= mps_state(H,d,4)
E_exact = expectation_value(ψ,H)

vals = Float64[]
n_exact= 3:8

for n in n_exact 
    problem = one_step_approx(h,n,SCS.Optimizer)
    push!(vals,problem.optval)
end

ΔELTI = real.(E_exact .- vals)

plt = Plots.plot(log2.(n_exact), log10.(ΔELTI), label="ELTI", xlabel="n", ylabel="log10(ΔE)", legend=:topleft)

D = 2
A = rand_unitary(d*D);
A = A[1:D,:];
A = reshape(A,D,d,D);
n_D2 = 5:40
vals = Float64[]
for n in n_D2
    problem = two_step_approx(h,D,n,A,SCS.Optimizer)
    push!(vals,problem.optval)
end

ΔErlxD2 = real.(E_exact .- vals)

plt = Plots.plot!(plt,log2.(n_D2), log10.(ΔErlxD2), label="E_relax_D2", xlabel="n", ylabel="log10(ΔE)", legend=:topleft)


D = 3
A = rand_unitary(d*D);
A = A[1:D,:];
A = reshape(A,D,d,D);
n_D3 = 5:60
vals = Float64[]
for n in n_D3
    problem = two_step_approx(h,D,n,A,SCS.Optimizer)
    push!(vals,problem.optval)
end

ΔErlxD2 = real.(E_exact .- vals)

plt = Plots.plot!(plt,log2.(n_D2), log10.(ΔErlxD2), label="E_relax_D2", xlabel="n", ylabel="log10(ΔE)", legend=:topleft)



δ = 1.0
H = heisenberg_XXZ(; Delta=δ , spin=1 // 2);
h = mat((kron(X, X) + kron(Y, Y) + δ*kron(Z, Z)) / 4);
ψ = mps_state(H, d, D)

A = rand_unitary(d*D);
A = A[1:D,:];
A = reshape(A,D,d,D);
typeof(A)

n =  20
res = main(h,D,n,A)

@show res.optval