using Test, RGRDMT, Yao
using MPSKit, MPSKitModels, TensorKit
using SCS, MosekTools, Dualization
using Random
using JuMP

function booda_dual()
    h = mat(-kron(Z, Z) / 4 - 1 / 4 * kron(X, I2))
    main_dual(h, 3:10, -1 / π, "data/etfi_dual.csv", "data/ntfi_dual.csv", SCS.Optimizer)
end

booda_dual() # this is currently the best option, could also try on Linux machinehttps://discourse.julialang.org/t/multi-core-parallel-support-for-jump-supported-solvers/112392/3 

function booda2()
    h = mat(-kron(Z, Z) / 4 - 1 / 4 * kron(X, I2))
    H = transverse_field_ising(; J=1 / 4)
    d = 2
    D = 2
    ψ_good = good_ground_state(H, 50)
    ψ = approx_ground_state(H, ψ_good, d, D)
    @test real(expectation_value(ψ, H)[]) ≈ -1 / π atol = 1e-3
    AL = ψ.AL[]

    k0 = Int(floor(2 * log(D) / log(d)) + 1)
    V0, L, R = CGmapping_from_AL(AL, k0, k0 + 2)

    main2(h, D, -1 / π, 8:12, V0, L, R, "data/etfi2.csv", "data/ntfi2.csv", MosekTools.Optimizer)
end

booda2()

function dooda_dual()
    h = mat(-kron(X, X) - kron(Y, Y) + kron(Z, Z)) / 4
    main_dual(h, 3:13, 1 / 4 - log(2), "data/exxx.csv", "data/nxxx.csv", MosekTools.Optimizer)
end

# dooda_dual()

function dooda2()
    h = mat(-kron(X, X) - kron(Y, Y) + kron(Z, Z)) / 4
    d = 2
    D = 4
    filename = "data/HeisenbergSpinHalf_SubLatticeRotation_Data.mat"
    ψ, upperBd = load_MPS(filename, D)
    AL = ψ.AL[]

    k0 = Int(floor(2 * log(D) / log(d)) + 1)
    V0, L, R = CGmapping_from_AL(AL, k0, k0 + 2)
    V0 = Complex.(V0)
    L = Complex.(L)
    R = Complex.(R)

    main2(h, D, 1 / 4 - log(2), 8:2:30, V0, L, R, "data/exxx$D.csv", "data/nxxx$D.csv", MosekTools.Optimizer)
end

dooda2()