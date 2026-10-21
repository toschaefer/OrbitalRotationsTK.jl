```@meta
CurrentModule = OrbitalRotationsTK
```

# Developer guide

This page is a map of the code: it names the parts, shows how they fit together and states
the conventions shared by all of them. The details (signatures, formulas, contracts) live in
the docstrings, collected in the [Code reference](@ref), and the code itself remains the
authoritative documentation. A change of the design updates this page in the same commit.
When in doubt, follow the nearest existing example.

## Nouns and verbs

OrbitalRotationsTK builds on [DFTK](https://docs.dftk.org)'s plane-wave infrastructure (e.g.
`PlaneWaveBasis`, FFTs, pseudopotentials), takes the overlap densities from
[PsiTK](https://github.com/toschaefer/PsiTK.jl) and leaves the unitary optimization to
[Lucon](https://github.com/toschaefer/Lucon.jl). It follows their pattern: *nouns* are plain
structs that carry either data or configuration, rather than both; *verbs* are functions
that take the data they need as arguments and the configuration as a dispatch argument.

| noun | role | lives in |
|---|---|---|
| functional, e.g. [`NPL`](@ref) | physics: *what* is optimized, i.e. the operators ``\sigma_F``, the weights `w` and the scalar function `h` | `src/functionals/<family>/<name>.jl` |
| family, e.g. [`JointDiagonalizationFunctional`](@ref) | physics: the common form of the loss of its functionals, and the numerics they share | `src/functionals/<family>/<family>.jl` |
| representation, e.g. [`OrbitalSubspace`](@ref) | numerics: *how* the loss and its gradient are computed | `src/common/representations.jl` |
| scalar function, e.g. [`Monomial`](@ref) | configuration: ``h`` with its derivative and Taylor degree | `src/common/scalar_functions.jl` |
| cache, e.g. [`JointDiagOrbitalSubspaceCache`](@ref) | data: what one combination of family and representation needs in every step, with its buffers | `src/functionals/<family>/<family>.jl` |
| [`RotationResult`](@ref) | data: the rotated orbitals, the optimal unitary and convergence information | `src/rotate.jl` |

The verbs:

| verb | transforms | lives in |
|---|---|---|
| [`rotate`](@ref) | orbitals + functional → [`RotationResult`](@ref) | `src/rotate.jl` |
| [`prepare_gradient`](@ref) | orbitals + functional + representation → cache | `src/functionals/<family>/<family>.jl` |
| [`gradient`](@ref) | cache + unitaries → Euclidean gradient and loss | `src/functionals/<family>/<family>.jl` |
| [`one_body_operators`](@ref) | functional + basis → operators ``\sigma_F`` | `src/functionals/<family>/<name>.jl` |
| [`maximize`](@ref), [`max_taylor_degree`](@ref) | functional → settings of the optimizer | `src/functionals/<family>/<family>.jl` |
| [`derivative`](@ref), [`taylor_degree`](@ref) | scalar function → ``h'(x)``, degree | `src/common/scalar_functions.jl` |

A typical call looks like this:

```julia
ψ   = DFTK.select_occupied_orbitals(scfres.basis, scfres.ψ, scfres.occupation).ψ  # DFTK
res = rotate(scfres.basis, ψ, NPL(w=[1.0, 0.5, 0.5]); representation=OrbitalSubspace())
```

[`rotate`](@ref) splits physics from numerics, so that functionals and representations
combine freely. It

1. checks its arguments and, unless `U0` is given, draws random initial unitaries;
2. calls [`prepare_gradient`](@ref)`(functional, representation, basis, ψ)` once, which
   dispatches on the family of the functional and on the representation, and asks the
   functional only for its operators via [`one_body_operators`](@ref);
3. passes [`maximize`](@ref) and [`max_taylor_degree`](@ref) of the functional to
   `Lucon.optimize`, which calls [`gradient`](@ref)`(prep, U, calc_loss) -> (Γ, L)` in
   every step;
4. returns a [`RotationResult`](@ref) with the rotated orbitals `ψ[k] * U[k]`.

`src/OrbitalRotationsTK.jl` includes the files in dependency order and declares the exports
(user-facing nouns and verbs) and `public` names (extension points) next to the `include`
of the file defining them. Every export needs a docstring (`docs/make.jl` builds with
`checkdocs = :exports`).

## Conventions

- **k-points:** per-k-point data are `Vector`s over k-points in the order of
  `basis.kpoints` (spin included), e.g. `ψ[k]` is an ``N_\text{pw} \times N`` matrix of
  orbital coefficients and `U[k]`, `Γ[k]` are ``N \times N``. Inside a cache, the orbitals
  of all k-points share one combined index ``(m, k)`` with the orbital running fastest.
  Lucon optimizes a single unitary, so [`rotate`](@ref) is Γ-point only so far.
- **Gradient:** [`gradient`](@ref) returns the Euclidean gradient
  ``\Gamma_{pq} = \partial L / \partial \overline{U_{pq}}``; Lucon turns it into the
  Riemannian one. The loss has the units stated in the docstring of each functional, which
  set the scale of `tol` in [`rotate`](@ref).
- **Preparation and steps:** [`prepare_gradient`](@ref) does all work that does not depend
  on `U` and allocates every buffer; [`gradient`](@ref) runs in every step and works in
  these buffers. The returned `Γ` may be a buffer that the next call overwrites.
- **Devices:** buffers are allocated with `similar` from the data they belong to, and
  scalar functions are applied by broadcasting, so a scalar function must be usable inside
  GPU kernels (see [`CustomFunction`](@ref)). The code is CPU-only so far.
- **Errors and assertions:** check a requirement in the routine whose algorithm needs
  it, as early as the user's choice reaches it (a public function, a constructor, or a
  validation helper called first). Throw a typed exception whose message says how to fix
  it: `DimensionMismatch` for inconsistent sizes, `ArgumentError` for any other invalid
  or unsupported input. Keep all "not supported (yet)" checks of a feature together in
  that one place. Low-level helpers carry no such policy and only `@assert` the
  preconditions their callers guarantee; a failing assertion is a bug in
  OrbitalRotationsTK, never a user error, and code must stay correct if asserts are not
  evaluated. E.g. [`rotate`](@ref) rejects bases beyond the Γ point, which
  [`prepare_gradient`](@ref) of [`JointDiagonalizationFunctional`](@ref) then only
  asserts, and the fallback of [`prepare_gradient`](@ref) rejects a combination of
  functional and representation that has no method.

## Extension points

Extending the code usually means adding one noun and its methods of the extension points.
The docstring of each extension point states its contract and what a new implementation
provides.

### Functionals and families

A family fixes the form of the loss and its gradient, and implements the numerics for each
representation it supports; its functionals only provide the physics. For the
[`JointDiagonalizationFunctional`](@ref) family, the physics are the one-body operators
``\sigma_F``.

- **A new functional of an existing family** is a configuration struct, a subtype of the
  family, with the fields and methods its family requires: for
  [`JointDiagonalizationFunctional`](@ref) the fields `h` and `w` and methods of
  [`one_body_operators`](@ref), in `src/functionals/joint_diagonalization/<name>.jl`. Its
  docstring states the units of the loss.
- **A new family** is a folder in `src/functionals/` with an abstract type, whose docstring
  states the form of the loss and its gradient and what a subtype provides, plus methods of
  [`maximize`](@ref) and [`max_taylor_degree`](@ref) for it. For each representation it
  supports, it adds a cache struct and methods of [`prepare_gradient`](@ref) and
  [`gradient`](@ref).

### Representations

**A new representation** is a configuration struct in `src/common/representations.jl`, plus
a cache struct and methods of [`prepare_gradient`](@ref) and [`gradient`](@ref) in each
family that supports it.

### Scalar functions

**A new scalar function** is a callable struct with methods of [`derivative`](@ref) and
[`taylor_degree`](@ref) in `src/common/scalar_functions.jl`.

## Computational cost

The cost of [`rotate`](@ref) is set by the combination of family and representation, not
by either alone, so it is documented per combination on the page
[Functionals and representations](@ref) instead of in the docstrings. A new family or
representation adds its row or column to the table "Available combinations" there and, for
each new combination, a row to the table "Cost and memory" and a paragraph to "Notes on the
combinations". A change of the scaling updates these entries. The columns are

- **preparation:** the leading terms of [`prepare_gradient`](@ref) in O-notation and the
  step they come from; it runs once per call of [`rotate`](@ref);
- **per gradient evaluation:** the leading terms of one call of [`gradient`](@ref); the
  optimizer calls it a few times per iteration;
- **memory:** the peak memory in complex numbers (16 bytes each), during the preparation
  and afterwards.

Counts use the following symbols. In code, the name is the symbol without the underscore of
its subscript (``N_G`` ↔ `NG`); counts without a symbol are snake_case.

| docs | code | meaning |
|---|---|---|
| ``N`` | `N` | orbitals per k-point, i.e. columns of `ψ[k]` |
| ``N_k`` | `Nk` | k-points (spin included) |
| ``N_F`` | `NF` | one-body operators of the functional with nonzero weight, e.g. atoms in [`NPL`](@ref) |
| ``N_\text{pw}`` | `Npw` | plane waves of the orbital basis |
| ``N_G`` | `NG` | plane waves of the overlap densities, set by `Ecut_ratio` |
| ``N_r`` | `Nr` | points of the FFT grid |
