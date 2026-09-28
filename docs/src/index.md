# OrbitalRotationsTK.jl

Optimal unitary rotations of a set of [DFTK.jl](https://dftk.org) orbitals with respect to a
user-chosen loss functional, for instance to localize them. Only the given orbitals are mixed;
the space they span is unchanged.

```julia
using DFTK, OrbitalRotationsTK

scfres = self_consistent_field(basis)
ψk     = scfres.ψ[1][:, 1:4]                  # the orbitals to rotate
res    = rotate(scfres.basis, ψk, NPL())
```

All public functions and types are listed in the [Code reference](@ref).
