# Contributing

Contributions via pull request are welcome. Please:

- Add tests for new functionality as `@testitem`s with a tag (see `test/runtests.jl`).
- Add a docstring to every exported name; the documentation build fails otherwise.
- Follow the existing code style ([Blue Style](https://github.com/invenia/BlueStyle),
  lines of about 92 characters).
- Keep pull requests focused on a single change or feature.

Questions? Open an issue.

## Code design

| concept | answers | lives in |
|---|---|---|
| **Functional**, e.g. `NPL` | *what* is optimized: operators σ_F, weights `w`, scalar function `h` | `src/functionals/<family>/<name>.jl` |
| **Family**, e.g. `OneBodyFunctional` | the numerics shared by all its functionals | `src/functionals/<family>/<family>.jl` |
| **Representation**: `OrbitalSubspace`, `FourierSpace`, `RealSpace` | *how* loss and gradient are computed | `src/common/representations.jl` |

`rotate(basis, ψ, NPL())` runs:

1. `prep = prepare_gradient(functional, representation, basis, ψ)`, once. It dispatches on
   (family, representation) and fetches the functional's operators via `one_body_operators`.
2. `gradient(prep, U, calc_loss) -> (Γ, L)`, called by Lucon in every step.

To add

- **a functional to an existing family:** one file with a struct (fields `h`, `w`) and
  `one_body_operators` for `FourierSpace` and `RealSpace`;
- **a representation:** a struct in `representations.jl`, plus `prepare_gradient`/`gradient`
  in each family that supports it;
- **a family:** a folder with an abstract type and its own `prepare_gradient`/`gradient`.
