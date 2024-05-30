using TensorKit,  MPSKit, MPSKitModels
using Yao
using Convex, MosekTools, SCS

function mps_state(d, D)
    H = heisenberg_XXX(; spin=1 // 2)
    random_data = TensorMap(rand, ComplexF64, ℂ^D * ℂ^d, ℂ^D)
    state = InfiniteMPS([random_data])
    groundstate, cache, delta = find_groundstate(state, H, GradientGrassmann(; maxiter=50))
    return groundstate
end

d = 2
D = 3 
ψ = mps_state(d, D)
A = ψ.AL[]

# A = ψ.AL[] # index D_in * d, D_out
# @tensor res[i, j] := A[k, l, j] * adjoint(A)[i, k, l]
# @assert permute(res,(1,),(2,)) ≈ id(domain(A)) 

# B = ψ.AR[] # index D_out , d * D_in 
# B = permute(B,(1,),(3,2))
# @tensor res[a,b] := B[a,i,j] * B'[i,j,b]
# @assert permute(res,(1,),(2,))  ≈ id(codomain(B))

iremain = Diagonal(ones(ComplexF64, D))

A
A[]
A' # this is valid input in (12) format
@tensor W2[i,k,j,l] := A'[i,a,k] * A'[a,j,l]
W2 = permute(W2,(1,3),(2,4))

R = A'
L = permute(A',(2,),(1,3)) 
@tensor W3[i,j,l,m,n] := A'[i,a,l] * A'[a,b,m] * A'[b,j,n]
W3 = permute(W3,(1,2),(3,4,5))

reshape(R[],d*D,D)
W2.data
kron(reshape(R[],D,d*D),iremain) 
# test eqn(11)
W2R = kron(reshape(R[],D,d*D),iremain) * kron(W2.data,mat(I2))
W2L = kron(iremain,reshape(L[],D,d*D)) * kron(mat(I2),W2.data)
@assert reshape(W3[],D^2,d^3) ≈ W2R
@assert reshape(W3[],D^2,d^3) ≈ W2L
reshape(W3[],D^2,d^3) - W2R
reshape(W3[],D^2,d^3) - W2L
W2R - W2L

@assert W2R ≈ W2L

function test_tl(W,A)


end


function main(D, n, A)
    # D = 2 # bond dimension (for coarse graining the states)
    # n = 10 # maximum site of reduced density matrix
    d = 2 # spin physical dimension
    # 1d translational invariant hamiltonian local term
    h = mat((kron(X, X) + kron(Y, Y) + kron(Z, Z)) / 4)

    ρ3 = ComplexVariable(d^3, d^3)
    ψ = mps_state(d, D)
    A = ψ.AL[]

    @tensor W2[j, l, i, k] := A[m,k,l] * adjoint(A)[j,m,i]

    # W2 = reshape(W2, D^2, d^2)
    W2 = reshape(permute(W2, (1,2),(3,4))[],D^2,d^2)

    R2 = reshape(permute(A,(3,),(1,2))[], D, d*D)
    L2 = reshape(A[],d*D,D)'

    iremain = Diagonal(ones(ComplexF64, d * D))
    i2mat = mat(I2)

    ωs = [ComplexVariable(d^2 * D^2, d^2 * D^2) for _ in 1:n-3]

    constraints = [
        tr(ρ3) == 1.0,
        ρ3 ⪰ 0,
        # [tr(ω) = 1.0 for ω in ωs]..., # automatically satisfied
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
res = main(4, 5)
@show res.optval