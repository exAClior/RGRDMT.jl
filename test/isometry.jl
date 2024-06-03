using RGRDMT
using Test
using MPSKit, MPSKitModels
using LinearAlgebra
using TensorKit

@testset "Isometry mapping" begin
    d = 2 # spin 1/2 system

    J = 1.0
    g = 1 / 2.0 # critical TFI
    H = transverse_field_ising(; J=J, g=g)
    k0 = rand(2:8)
    Dmax = floor(Int, sqrt(d^k0 - 1))
    D = rand(2:Dmax)
    ψ, _, _ = approx_ground_state(H, d, D)

    A = ψ.AL[]

    @show D, k0
    # why do we have such long compilation time?
    V0, L, R = CGmapping_from_AL(A, k0)
    i2mat = diagm(ones(eltype(A), d))
    @test V0 * V0' ≈ diagm(ones(eltype(A), D^2))
    @test L * L' ≈ diagm(ones(eltype(A), D^2))
    @test R * R' ≈ diagm(ones(eltype(A), D^2))
    @test L * kron(i2mat, V0) ≈ R * kron(V0, i2mat)
end

@testset "upper bound ansatz" begin
    d = 2
    H = heisenberg_XXX(; spin=1 // 2)

    ψ_good = good_ground_state(H,100)
    real(expectation_value(ψ_good, H))

    D = 2
    ψ_approx = approx_ground_state(H,ψ_good, d, D)
    E_approx = real(expectation_value(ψ_approx, H))
    # E_exact = 0.25 - log(2)
end
