using Test, RGRDMT, Yao
using MPSKit, MPSKitModels, TensorKit
using SCS, MosekTools, Dualization
using Random
using JuMP

function booda_dual()
    h = mat(-kron(Z, Z) / 4 - 1 / 8 * kron(X, I2) - 1 / 8 * kron(I2, X))
    main_dual(h, 3:8, -1 / π, "data/etfi_dual.csv", "data/ntfi_dual.csv", SCS.Optimizer)
end

# booda_dual() # this is currently the best option, could also try on Linux machinehttps://discourse.julialang.org/t/multi-core-parallel-support-for-jump-supported-solvers/112392/3 

function booda2(D)
    h = mat(kron(Z, Z) / 4 + 1 / 8 * kron(X, I2) + 1 / 8 * kron(I2, X))
    H = transverse_field_ising(; J=1 / 4)
    d = 2
    ψ = good_ground_state(H, D)
    AL = ψ.AL[]

    k0 = Int(floor(2 * log(D) / log(d)) + 1)
    V0, L, R = CGmapping_from_AL(AL, k0, k0 + 2)

    main2(h, D, -1 / π, vcat(k0+2:2:k0+10, k0+12:10:70), V0, L, R, "data/etfi$D.csv", "data/ntfi$D.csv", MosekTools.Optimizer)
end

booda2(2)
booda2(3)
booda2(4)
booda2(5)
booda2(6)

function dooda_dual()
    h = mat(-kron(X, X) - kron(Y, Y) + kron(Z, Z)) / 4
    main_dual(h, 3:12, 1 / 4 - log(2), "data/exxx.csv", "data/nxxx.csv", SCS.Optimizer)
end

# dooda_dual()

function dooda2(D)
    h = mat(-kron(X, X) - kron(Y, Y) + kron(Z, Z)) / 4
    d = 2
    filename = "data/HeisenbergSpinHalf_SubLatticeRotation_Data.mat"
    ψ, upperBd = load_MPS(filename, D)
    AL = ψ.AL[]

    k0 = Int(floor(2 * log(D) / log(d)) + 1)
    V0, L, R = CGmapping_from_AL(AL, k0, k0 + 2)
    V0 = Complex.(V0)
    L = Complex.(L)
    R = Complex.(R)

    main2(h, D, 1 / 4 - log(2), 8:2:64, V0, L, R, "data/exxx$D.csv", "data/nxxx$D.csv", MosekTools.Optimizer)
end

# dooda2(2)
# dooda2(3)
# dooda2(4)
# dooda2(5)
# dooda2(6)
# dooda2(7)