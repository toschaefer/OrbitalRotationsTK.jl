# OrbitalRotationsTK.jl

[![Build Status](https://github.com/toschaefer/OrbitalRotationsTK.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/toschaefer/OrbitalRotationsTK.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://toschaefer.github.io/OrbitalRotationsTK.jl/dev/)

Optimal unitary rotations of a set of [DFTK.jl](https://dftk.org) orbitals with respect to a
user-chosen loss functional, for instance to localize them. Only the given orbitals are mixed
(occupied, virtual, or both); unlike an SCF step, the rotation leaves the space they span
unchanged. The unitary optimization is done by [Lucon.jl](https://github.com/toschaefer/Lucon.jl).

> **Status:** early development. The interface below is in place; the numerics are not
> implemented yet.

## Features

- Functionals: Nuclear Potential Localization, `NPL(; h, w)`, with a scalar function `h`
  (e.g. `Monomial(2)`) and optional weights `w` per atom.
- Representations of the loss and its gradient: `OrbitalSubspace()` (default).
- Γ-point only.

## Installation

OrbitalRotationsTK currently depends on development branches of DFTK.jl and PsiTK.jl,
declared in `[sources]` of `Project.toml`. Install it with `develop`, which honors them:

```julia
using Pkg
Pkg.develop(url="https://github.com/toschaefer/OrbitalRotationsTK.jl")
```

Once the required changes are released in DFTK.jl and PsiTK.jl, a plain `Pkg.add` will do.

## Usage

```julia
using DFTK, OrbitalRotationsTK

scfres = self_consistent_field(basis)
ψ      = scfres.ψ[1][:, 1:4]                  # the orbitals to rotate
res    = rotate(scfres.basis, ψ, NPL())

res.U   # optimal unitary
res.ψ   # rotated orbitals
```

## To do

- [ ] Implement `rotate` and the `OrbitalSubspace` evaluation for `NPL`
- [ ] `FourierSpace` and `RealSpace` representations
- [ ] More functionals: Foster–Boys, Pipek–Mezey, Edmiston–Ruedenberg, intrinsic bond orbitals
- [ ] k-points
- [ ] Switch from the development branches to released versions of DFTK.jl and PsiTK.jl

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md), including a short overview of the code design.
