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
    D = rand(2:7)
    ψ_good = good_ground_state(H, 30)
    ψ = approx_ground_state(H, ψ_good, d, D)
    @test real(expectation_value(ψ, H)[]) ≈ -1 / π atol = 1e-3
    AL = ψ.AL[]

    V0, L, R = CGmapping_from_AL(AL, k0)

    main2(h, D, -1 / π, 5:12, V0, L, R, "data/etfi2.csv", "data/ntfi2.csv", dual_optimizer(MosekTools.Optimizer))
end

function dooda_dual()
    h = mat(-kron(X, X) - kron(Y, Y) + kron(Z, Z)) / 4
    main(h, 3:10, 1 / 4 - log(2), "data/exxx.csv", "data/nxxx.csv", SCS.Optimizer)
end

dooda_dual()

function dooda2()
    h = mat(-kron(X, X) - kron(Y, Y) + kron(Z, Z)) / 4
    d = 2
    filename = "data/HeisenbergSpinHalf_SubLatticeRotation_Data.mat"
    D = rand(2:7)
    ψ, upperBd = load_MPS(filename, D)
    AL = ψ.AL[]

    V0, L, R = CGmapping_from_AL(AL, k0)
    V0 = Complex.(V0)
    L = Complex.(L)
    R = Complex.(R)

    main2(h, D, 1 / 4 - log(2), 5:12, V0, L, R, "data/exxx2.csv", "data/nxxx2.csv", dual_optimizer(MosekTools.Optimizer))
end

dooda2()