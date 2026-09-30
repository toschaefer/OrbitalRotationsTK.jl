# OrbitalRotationsTK.jl

Optimal unitary rotations of a set of [DFTK.jl](https://dftk.org) orbitals with respect to a
user-chosen loss functional, for instance to localize them. Only the given orbitals are
mixed; the space they span is unchanged.

## Installation

OrbitalRotationsTK currently depends on development branches of DFTK.jl and PsiTK.jl, so
install it with `develop`:

```julia
using Pkg
Pkg.develop(url="https://github.com/toschaefer/OrbitalRotationsTK.jl")
```

## Example

Localize the occupied orbitals of a converged DFTK calculation at the Γ point:

```julia
using DFTK, OrbitalRotationsTK

ψ   = [ψk[:, occk .> 0] for (ψk, occk) in zip(scfres.ψ, scfres.occupation)]
res = rotate(scfres.basis, ψ, NPL())
res.ψ    # the localized orbitals
```

## Where to go next

- The [Tutorial](@ref "Tutorial: localized orbitals of water") works through a complete
  example, a water molecule, and checks the localized orbitals.
- The [Code reference](@ref) lists the public API and, for developers, the documented
  internals.
