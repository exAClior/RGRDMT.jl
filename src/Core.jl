using TensorKit,  MPSKit, MPSKitModels
using Yao
using Convex, MosekTools, SCS

function mps_state(H,d, D)
    random_data = TensorMap(rand, ComplexF64, ℂ^D * ℂ^d, ℂ^D)
    state = InfiniteMPS([random_data])
    groundstate, cache, delta = find_groundstate(state, H, GradientGrassmann(; maxiter=50))
    return groundstate
end

function main(h, D, n, A)
    d = 2 # spin physical dimension
    # 1d translational invariant hamiltonian local term

    ρ3 = ComplexVariable(d^3, d^3)

    @tensor W2[j, l, i, k] := A[j,i,a] * A[a,k,l]

    W2 = reshape(W2, D^2, d^2)

    R2 = reshape(A, D,D*d) 
    L2 = transpose(reshape(A,D*d,D))

    iremain = Diagonal(ones(ComplexF64,d* D))
    i2mat = mat(I2)

    ωs = [ComplexVariable(d^2 * D^2, d^2 * D^2) for _ in 1:n-3]

    constraints = [
        tr(ρ3) == 1.0,
        ρ3 ⪰ 0,
        [ω ⪰ 0 for ω in ωs]...,
        partialtrace(ρ3, 1, d * ones(Int64, 3)) == partialtrace(ρ3, 3, d * ones(Int64, 3)),
        kron(W2, i2mat) * ρ3 * kron(W2, i2mat)' == partialtrace(ωs[1], 1, [d, D^2, d]),
        kron(i2mat, W2) * ρ3 * kron(i2mat, W2)' == partialtrace(ωs[1], 3, [d, D^2, d])
    ]

    for ii in 2:n-3
        push!(constraints, kron(L2, iremain) * ωs[ii-1] * kron(L2, iremain)' == partialtrace(ωs[ii], 1, [d, D^2, d]))
        push!(constraints, kron(iremain, R2) * ωs[ii-1] * kron(iremain, R2)' == partialtrace(ωs[ii], 3, [d, D^2, d]))
    end

    problem = minimize(real(tr(kron(h, mat(I2)) * ρ3)), constraints)

    solve!(problem, Mosek.Optimizer)
    return problem
end


d = 2
D = 4 

δ = 0.5
H = heisenberg_XXZ(; Delta=δ , spin=1 // 2)
h = mat((kron(X, X) + kron(Y, Y) + δ*kron(Z, Z)) / 4)
ψ = mps_state(H, d, D)

A = rand_unitary(d*D)
A = A[1:D,:]
A = reshape(A,D,d,D)

n =  5 
res = main(h,D,n,A)

@show res.optval