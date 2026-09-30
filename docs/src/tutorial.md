# Tutorial: localized orbitals of water

This tutorial localizes the four occupied orbitals of a water molecule with the
[`NPL`](@ref) functional. It runs in well under a minute. Besides OrbitalRotationsTK,
it needs [PseudoPotentialData.jl](https://github.com/JuliaMolSim/PseudoPotentialData.jl)
for the pseudopotentials.

## Ground state

A water molecule in a cubic box, computed with DFTK at the Γ point only
(`kgrid=(1, 1, 1)`), as [`rotate`](@ref) currently requires:

```julia
using DFTK, OrbitalRotationsTK, PseudoPotentialData, LinearAlgebra

pseudopotentials = PseudoFamily("dojo.nc.sr.lda.v0_4_1.standard.upf")
O = ElementPsp(:O, pseudopotentials)
H = ElementPsp(:H, pseudopotentials)

a = 12.0                                              # cubic box (bohr)
lattice = a * I(3)
positions_cart = [[0.0, 0.0, 0.0], [1.43, 1.11, 0.0], [-1.43, 1.11, 0.0]]
positions = [r / a .+ 0.5 for r in positions_cart]    # reduced coordinates, box centre

model  = model_DFT(lattice, [O, H, H], positions; functionals=LDA())
basis  = PlaneWaveBasis(model; Ecut=20, kgrid=(1, 1, 1))
scfres = self_consistent_field(basis; tol=1e-8)
```

## Localization

DFTK computes a few unoccupied orbitals in addition to the occupied ones, so we first select
the occupied orbitals, one matrix per k-point. [`rotate`](@ref) then finds the unitary that
maximizes the NPL functional:

```julia
ψ   = DFTK.select_occupied_orbitals(scfres.basis, scfres.ψ, scfres.occupation).ψ
res = rotate(scfres.basis, ψ, NPL())
```

`res.ψ` holds the rotated orbitals, again one matrix per k-point, and `res.optimizer` the
result of the optimization, e.g. the optimal unitary `res.optimizer.U` and the final loss
`res.optimizer.loss`. For water, the rotated orbitals are the two O–H bonds and two lone
pairs on oxygen. The optimization starts from a random unitary; pass
`rng=Random.Xoshiro(1)` to [`rotate`](@ref) for identical results in every run.

## Options

- **The functional:** `NPL(h=Monomial(4))` changes the scalar function ``h``, and
  `NPL(w=[1.0, 0.5, 0.5])` weights the atoms, in the order of `basis.model.atoms`; a weight
  of zero excludes an atom. See [`NPL`](@ref).
- **The optimization:** `tol`, `maxiter`, `U0` and `callback`, e.g.
  `callback=OrbitalRotationsTK.Lucon.PrintTrace()` for a convergence trace. See
  [`rotate`](@ref).
