```@meta
CurrentModule = OrbitalRotationsTK
```

# Functionals and representations

[`rotate`](@ref) evaluates the loss and its gradient through a combination of two choices:

- the **functional**, e.g. [`NPL`](@ref): *what* is optimized. Functionals of the same form
  belong to one family, e.g. [`JointDiagonalizationFunctional`](@ref), and share the
  numerics.
- the **representation**, e.g. [`OrbitalSubspace`](@ref), passed as
  `rotate(...; representation)`: *how* the loss and its gradient are computed.

The cost and memory depend on the combination of family and representation. This page
lists which combinations exist and how they scale.

## Available combinations

| family (functionals)                                   | `OrbitalSubspace` | `FourierSpace`        | `RealSpace`           |
|:-------------------------------------------------------|:------------------|:----------------------|:----------------------|
| [`JointDiagonalizationFunctional`](@ref) ([`NPL`](@ref)) | available         | not yet implemented   | not yet implemented   |

## Cost and memory

Symbols: ``N`` orbitals, ``N_F`` operators of the functional (e.g. atoms with nonzero
weight in [`NPL`](@ref)), ``N_G`` plane waves of the overlap densities (see `Ecut_ratio`),
``N_r`` points of the FFT grid. The preparation runs once per call of [`rotate`](@ref); the
optimizer then evaluates the gradient a few times per iteration. All numbers are for the Γ
point, currently the only supported k-point.

| family                          | representation    | preparation                                            | per gradient evaluation   | memory                                                   |
|:--------------------------------|:------------------|:-------------------------------------------------------|:--------------------------|:---------------------------------------------------------|
| `JointDiagonalizationFunctional` | `OrbitalSubspace` | ``O(N^2 N_r \log N_r)`` FFTs, ``O(N^2 N_G N_F)`` contraction | ``O(N_F N^3)``           | ``N^2 N_G`` during the preparation, then ``N^2 N_F``      |

Memory is given in complex numbers (16 bytes each).

## Notes on the combinations

**`JointDiagonalizationFunctional` with `OrbitalSubspace`:** the matrix elements
``\langle \psi_m | \sigma_F | \psi_n \rangle`` are computed once from the overlap densities
of all orbital pairs. Afterwards, each gradient evaluation is independent of the plane-wave
basis and consists of ``N_F`` matrix products of size ``N \times N``. The practical limit is
the memory of the overlap densities during the preparation.
