# OrbitalRotationsTK.jl

Optimal unitary rotations of a set of [DFTK.jl](https://dftk.org) orbitals with respect to a
user-chosen loss functional, for instance to localize them. Only the given orbitals are mixed;
the space they span is unchanged.

```julia
using DFTK, OrbitalRotationsTK

scfres = self_consistent_field(basis)
ψ      = [ψk[:, 1:4] for ψk in scfres.ψ]      # the orbitals to rotate, per k-point
res    = rotate(scfres.basis, ψ, NPL())
```

All public functions and types are listed in the [Code reference](@ref).
