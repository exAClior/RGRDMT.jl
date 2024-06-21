using RGRDMT, Test
using TensorKit, MPSKit, MPSKitModels
using MAT

@testset "Hamiltonian" begin
    vars = matread("HeisenbergSpinHalf_SubLatticeRotation_Data.mat")
    h_matlab = vars["H"] # why is their hamiltonian equal to (ZZ - XX - YY) /4 
    Sx = TensorMap([0 1; 1 0] / 2, ℂ^2, ℂ^2)
    Sy = TensorMap([0 -im; im 0] / 2, ℂ^2, ℂ^2)
    Sz = TensorMap([1 0; 0 -1] / 2, ℂ^2, ℂ^2)
    h_tensorkit = (-Sx ⊗ Sx - Sy ⊗ Sy + Sz ⊗ Sz)
    @test h_matlab ≈ reshape(h_tensorkit[], 4, 4) # not sure why, their hamiltonian is actually (ZZ - XX - YY) /4
end

@testset "Load MPS" begin
    filename = "HeisenbergSpinHalf_SubLatticeRotation_Data.mat"
    ham_MPO = heisenberg_XXZ(; J=-1.0, Delta=-1.0, spin=1 // 2)
    for D in 2:8
        ψ, upperBd, E_exact = load_MPS(filename, D)
        AL = ψ.AL[]
        @tensor res[-1; -2] := AL[1 2 -2] * conj(AL[1 2 -1])
        @test res ≈ id(space(res, 1))

        AR = ψ.AR[]
        @tensor res[-1; -2] := AR[-1 2 1] * conj(AR[-2 2 1])
        @test res ≈ id(space(res, 1))

        AC = ψ.AC[]
        C = ψ.CR[]
        @test AL * C ≈ AC

        @test upperBd ≈ expectation_value(ψ, ham_MPO)[1] atol = 1e-14
    end
end