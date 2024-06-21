using TensorKit, MPSKit, MPSKitModels
using Yao

# d = 2
# D = 137
# A0 = TensorMap(reshape(rand_unitary(d * D)[:, 1:D], D, d, D), ℂ^D * ℂ^d, ℂ^D)
# ψ0 = InfiniteMPS([A0])
# H0 = heisenberg_XXX(; J=4.0, spin=1 // 2)
# ψ_gs, _ = find_groundstate(ψ0, H0, VUMPS());
# expectation_value(ψ_gs, H0)

# ψ0 = InfiniteMPS([ℂ^2], [ℂ^20])
# H0 = transverse_field_ising(; J=1 / 4)
# ψ_gs, _ = find_groundstate(ψ0, H0, VUMPS());
# expectation_value(ψ_gs, H0)

# tfi analytic is 1/π%