```@meta
CurrentModule = RGRDMT
```

# Renormalization Group Reduced Density Matrix 


RGRDMT.jl is a julia library for finding the lower bound to the ground state energy of a translationally invariant Hamiltonian. The library is based on the paper [^1].

[^1]:[kull2024lower](@cite) 

This library uses Convex Optimization[^2] to find the lower bound to the ground state energy by optimizing a reduced density matrix. Concepts of renormalization group[^3] were employed to reduce the dimension of the problem.

[^2]:[boyd2004convex](@cite)
[^3]:[xiang2023density](@cite)

## Quick Start

Currently, we have included tested optimizers for Convex Optimization and the package is plug-and-play. You will only need to install the package by running the following code.

```julia
using Pkg
Pkg.add("RGRDMT")
```

For detailed usage, please refer to the test case in `test/example.jl`

Documentation for [RGRDMT](https://github.com/exAClior/RGRDMT.jl).


