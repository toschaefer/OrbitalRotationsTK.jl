# Contributing

Contributions via pull request are welcome. Please:

- Add tests for new functionality as `@testitem`s with a tag (see `test/runtests.jl`).
- Add a docstring to every exported name; the documentation build fails otherwise.
- Follow the existing code style ([Blue Style](https://github.com/invenia/BlueStyle),
  lines of about 92 characters).
- Keep pull requests focused on a single change or feature.

Questions? Open an issue.

## Code design

This section is a short orientation, not a reference: it shows the pattern the code follows
through a few examples, so that new code can follow it too. The docstrings and the code
itself remain the authoritative documentation.

OrbitalRotationsTK builds on DFTK's plane-wave infrastructure (e.g. `PlaneWaveBasis`, FFTs,
pseudopotentials), takes overlap densities from PsiTK and leaves the unitary optimization to
Lucon. It follows their pattern: *nouns* are plain structs that carry either data or
configuration, rather than both; *verbs* are functions that take the data as arguments and
the configuration as a dispatch argument. Some examples:

| noun | role | lives in |
|---|---|---|
| functional, e.g. `NPL` | physics: *what* is optimized, i.e. the operators σ_F, weights `w` and scalar function `h` | `src/functionals/<family>/<name>.jl` |
| family, e.g. `OneBodyFunctional` | physics: the common form of its functionals, and the numerics they share | `src/functionals/<family>/<family>.jl` |
| representation, e.g. `OrbitalSubspace` | numerics: *how* the loss and its gradient are computed | `src/common/representations.jl` |
| scalar function, e.g. `Monomial` | configuration: `h` with its derivative and Taylor degree | `src/common/scalar_functions.jl` |
| `RotationResult` | data: the optimal unitary, the rotated orbitals, convergence information | `src/rotate.jl` |

A typical call looks like this:

```julia
ψ   = scfres.ψ[1][:, 1:4]                                    # DFTK enters here
res = rotate(scfres.basis, ψ, NPL(w=[1.0, 0.5, 0.5]); representation=OrbitalSubspace())
```

`rotate` illustrates the split between physics and numerics:
`prepare_gradient(functional, representation, basis, ψ)` dispatches on the family of the
functional and on the representation, and asks the functional only for its operators via
`one_body_operators`; Lucon then calls `gradient(prep, U, calc_loss) -> (Γ, L)` in every
step. Functionals and representations therefore combine freely.

Extending the code usually means adding one such noun and its methods, for instance

- **a functional of an existing family:** a configuration struct with the fields `h` and
  `w`, plus `one_body_operators` methods for `FourierSpace` and `RealSpace`, in
  `src/functionals/one_body/<name>.jl`;
- **a representation:** a configuration struct in `src/common/representations.jl`, plus
  `prepare_gradient` and `gradient` methods in each family that supports it;
- **a family with a different form of the loss:** a folder in `src/functionals/` with an
  abstract type and its own `prepare_gradient` and `gradient` methods;
- **a scalar function:** a struct with `derivative` and `taylor_degree` methods in
  `src/common/scalar_functions.jl`.

When in doubt, follow the nearest existing example.
