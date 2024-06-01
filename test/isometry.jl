using RGRDMT
using Test
using MPSKit, MPSKitModels
using LinearAlgebra
using TensorKit

@testset "Isometry mapping" begin
    J = 1.0
    g = 1 / 2.0 # critical TFI
    H = transverse_field_ising(; J=J, g=g)
    D = 3
    d = 2
    ψ, _, _ = mps_state(H, d, D);

    A = ψ.AL[]

    V0, L, R = CGmapping_from_AL(A, 4)
    i2mat = diagm(ones(eltype(A), d))
    @test V0 * V0' ≈ diagm(ones(eltype(A),D^2))
    @test L * L' ≈ diagm(ones(eltype(A),D^2))
    @test R * R' ≈ diagm(ones(eltype(A),D^2))

    @test L * kron(i2mat,V0) ≈ R * kron(V0,i2mat) 
end