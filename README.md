# OrbitalRotationsTK.jl

[![Build Status](https://github.com/toschaefer/OrbitalRotationsTK.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/toschaefer/OrbitalRotationsTK.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://toschaefer.github.io/OrbitalRotationsTK.jl/dev/)

Optimal unitary rotations of a set of [DFTK.jl](https://dftk.org) orbitals with respect to a
user-chosen loss functional, for instance to localize them. Only the given orbitals are mixed
(occupied, virtual, or both); unlike an SCF step, the rotation leaves the space they span
unchanged. The unitary optimization is done by [Lucon.jl](https://github.com/toschaefer/Lucon.jl).

> **Status:** early development. The `NPL` functional works at the Γ point; the interface
> may still change.

## Features

- Functionals: Nuclear Potential Localization, `NPL(; h, w)`, with a scalar function `h`
  (e.g. `Monomial(2)`) and optional weights `w` per atom.
- Representations of the loss and its gradient: `OrbitalSubspace()` (default). The
  [available combinations and their cost](https://toschaefer.github.io/OrbitalRotationsTK.jl/dev/functionals_and_representations/)
  are listed in the documentation.
- Starts from a random unitary, so that symmetric structures do not get stuck at a saddle
  point; reproducible with a seeded `rng`.
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

Localize the occupied orbitals of a converged DFTK calculation at the Γ point:

```julia
using DFTK, OrbitalRotationsTK

ψ   = DFTK.select_occupied_orbitals(scfres.basis, scfres.ψ, scfres.occupation).ψ
res = rotate(scfres.basis, ψ, NPL())

res.ψ             # the localized orbitals, one matrix per k-point
res.optimizer.U   # the optimal unitary
```

The [documentation](https://toschaefer.github.io/OrbitalRotationsTK.jl/dev/) has a
[tutorial](https://toschaefer.github.io/OrbitalRotationsTK.jl/dev/tutorial/) with a complete
example and the reference of all functions.

## To do

- [x] Implement `rotate` and the `OrbitalSubspace` evaluation for `NPL`
- [ ] Tests for `gradient` and `rotate`
- [ ] `FourierSpace` and `RealSpace` representations
- [ ] More functionals: Foster–Boys, Pipek–Mezey, Edmiston–Ruedenberg, intrinsic bond orbitals
- [ ] k-points
- [ ] Switch from the development branches to released versions of DFTK.jl and PsiTK.jl

## Contributing and Support

We welcome contributions from the scientific community! 

- If you encounter a bug, have a feature request, or need help, please open an [issue](https://github.com/toschaefer/OrbitalRotationsTK.jl/issues).
- If you'd like to contribute code, please submit a Pull Request. We recommend opening an issue first to discuss your planned changes. [CONTRIBUTING.md](CONTRIBUTING.md) explains the code design and how to extend it.
