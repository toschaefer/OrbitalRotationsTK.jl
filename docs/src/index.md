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

ψ   = DFTK.select_occupied_orbitals(scfres.basis, scfres.ψ, scfres.occupation).ψ
res = rotate(scfres.basis, ψ, NPL())
res.ψ    # the localized orbitals
```

## Where to go next

- The [Tutorial](@ref "Tutorial: localized orbitals of water") works through a complete
  example, a water molecule.
- [Functionals and representations](@ref) shows which combinations are available and how
  their cost and memory scale with the system size.
- The [Code reference](@ref) lists the public API and, for developers, the documented
  internals.
