using Test, RGRDMT, Yao
using MPSKit, MPSKitModels, TensorKit
using SCS, MosekTools, Dualization
using Random
using JuMP

# @testset "Performance" begin
#     using BenchmarkTools
#     n = 6
#     optimizer = MosekTools.Optimizer
#     h = mat(-kron(Z, Z) / 4 - 1 / 4 * kron(X, I2))
#     problem = one_step_approx(h, n, dual_optimizer(optimizer)) # 7.74s 43.338MiB
#     problem = one_step_approx(h, n, optimizer) # 20s didn't terminate 66.720MiB
#     model = one_step_approx_jp(h, n, optimizer) # 13.53s 115212272 bytes
#     model = one_step_approx_jp(h, n, dual_optimizer(optimizer)) # 10.62s
# end

function booda_dual()
    h = mat(-kron(Z, Z) / 4 - 1 / 4 * kron(X, I2))
    main_dual(h, 3:10, -1 / π, "data/etfi_dual.csv", "data/ntfi_dual.csv", MosekTools.Optimizer)
end

booda_dual()

function booda()
    h = mat(-kron(Z, Z) / 4 - 1 / 4 * kron(X, I2))
    main(h, 3:10, -1 / π, "data/etfi.csv", "data/ntfi.csv", dual_optimizer(MosekTools.Optimizer))
end

# booda()

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

# booda2()

function dooda()
    h = mat(-kron(X, X) - kron(Y, Y) + kron(Z, Z)) / 4
    main(h, 11:13, 1 / 4 - log(2), "data/exxx.csv", "data/nxxx.csv", dual_optimizer(MosekTools.Optimizer))
end

dooda()

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