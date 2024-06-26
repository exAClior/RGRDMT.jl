using RGRDMT, Test
using MPSKit, MPSKitModels, TensorKit
using LinearAlgebra
using SparseArrays

@testset "Isometry mapping tfi" begin
    d = 2 # spin 1/2 system

    J = 1 / 4
    H = transverse_field_ising(; J=J)
    D = rand(2:7)
    ψ = good_ground_state(H, D)

    AL = ψ.AL[]

    k0 = Int(floor(2 * log(D) / log(d)) + 1)
    n = k0 + 2
    V0, L, R = CGmapping_from_AL(AL, k0, n)
    sp_eyed = sparse(I, d, d)
    @test R * R' ≈ diagm(ones(eltype(AL), D^2))
    @test L * kron(sp_eyed, V0) ≈ R * kron(V0, sp_eyed)
end


@testset "Isometry mapping xxx" begin
    d = 2 # spin 1/2 system

    filename = "data/HeisenbergSpinHalf_SubLatticeRotation_Data.mat"
    filename = joinpath(dirname(dirname(@__FILE__)), filename)
    H = heisenberg_XXZ(; J=-1.0, Delta=-1.0, spin=1 // 2)

    D = rand(2:7)
    ψ, upperBd = load_MPS(filename, D)
    AL = ψ.AL[]

    k0 = Int(floor(2 * log(D) / log(d)) + 1)
    n = k0 + 2

    V0, L, R = CGmapping_from_AL(AL, k0, n)
    sp_eyed = sparse(I, d, d)
    @test R * R' ≈ diagm(ones(eltype(AL), D^2)) # since we are using AL
    @test L * kron(sp_eyed, V0) ≈ R * kron(V0, sp_eyed)
end

# test is too long, only run it locally
# @testset "TFI Model Ansatz" begin
#     d = 2
#     D = 20
#     H = transverse_field_ising(; J=1 / 4)

#     ψ_good = good_ground_state(H, D)

#     E_exact = -1 / π # there was a typo in pfeuty 1970, use lieb 1961
#     @test real(expectation_value(ψ_good, H)[]) ≈ E_exact atol = 1e-6

#     D = 2
#     ψ_approx = approx_ground_state(H, ψ_good, d, D)
#     @test real(expectation_value(ψ_approx, H)[]) ≈ E_exact atol = 1e-3
# end

# test takes too long, only run it locally
# @testset "XXZ Model Ansatz" begin
#     d = 2
#     D = 40 # change this to 137, you should see fast dropping in gradient
#     H = heisenberg_XXZ(; J=-1.0, Delta=-1.0, spin=1 // 2)

#     ψ_good = good_ground_state(H, D)

#     E_exact = 0.25 - log(2)
#     @test real(expectation_value(ψ_good, H)[]) ≈ E_exact atol = 1e-8

#     # approximation is no good don't use
#     D = 2
#     ψ_approx = approx_ground_state(H, ψ_good, d, D)
#     _, E_inpaper = load_MPS("data/HeisenbergSpinHalf_SubLatticeRotation_Data.mat", D)
#     @test real(expectation_value(ψ_approx, H)[]) ≈ E_inpaper
# end
