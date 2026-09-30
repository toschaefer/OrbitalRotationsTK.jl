# Tutorial: localized orbitals of water

This tutorial localizes the four occupied orbitals of a water molecule with the
[`NPL`](@ref) functional and checks that the result consists of two O–H bonds and two lone
pairs. It runs in well under a minute. Besides OrbitalRotationsTK, it needs
[PseudoPotentialData.jl](https://github.com/JuliaMolSim/PseudoPotentialData.jl) for the
pseudopotentials.

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

DFTK computes a few unoccupied orbitals in addition to the occupied ones, so we select the
occupied orbitals by their occupation, one matrix per k-point. [`rotate`](@ref) then finds
the unitary that maximizes the NPL functional:

```julia
ψ   = [ψk[:, occk .> 0] for (ψk, occk) in zip(scfres.ψ, scfres.occupation)]
res = rotate(scfres.basis, ψ, NPL())
```

`res.ψ` holds the rotated orbitals, again one matrix per k-point, and `res.optimizer` the
result of the optimization, e.g. the optimal unitary `res.optimizer.U` and the final loss
`res.optimizer.loss`. The optimization starts from a random unitary; pass
`rng=Random.Xoshiro(1)` to [`rotate`](@ref) for identical results in every run.

## Checking the result

The centre of each rotated orbital, and its distance to the three atoms:

```julia
r = vec(r_vectors_cart(basis))
atoms_cart = [lattice * p for p in positions]
for (i, ϕ) in enumerate(eachcol(res.ψ[1]))
    ρ = vec(abs2.(ifft(basis, basis.kpoints[1], ϕ)))
    centre = sum(ρ .* r) / sum(ρ)
    distances = [norm(centre - R) for R in atoms_cart]
    println("orbital $i: distances to O, H, H = ", round.(distances; digits=2))
end
```

```
orbital 1: distances to O, H, H = [0.09, 1.75, 1.75]
orbital 2: distances to O, H, H = [0.8, 1.02, 2.08]
orbital 3: distances to O, H, H = [0.8, 2.08, 1.02]
orbital 4: distances to O, H, H = [0.45, 2.12, 2.12]
```

The order of the orbitals may differ between runs.

- **Orbitals 2 and 3 are the two O–H bonds:** their centres lie 0.80 bohr from O and
  1.02 bohr from one H, which adds up to the O–H bond length of 1.81 bohr.
- **Orbitals 1 and 4 are the lone pairs on oxygen:** one lies in the plane of the molecule,
  the other perpendicular to it. They are not equivalent, as known from atom-based
  localization such as Pipek–Mezey.

The canonical orbitals `ψ`, in contrast, are delocalized: each of them has its centre on
the symmetry axis of the molecule, equally far from both hydrogens.

## Options

- **The functional:** `NPL(h=Monomial(4))` changes the scalar function ``h``, and
  `NPL(w=[1.0, 0.5, 0.5])` weights the atoms, in the order of `basis.model.atoms`; a weight
  of zero excludes an atom. See [`NPL`](@ref).
- **The optimization:** `tol`, `maxiter`, `U0` and `callback`, e.g.
  `callback=OrbitalRotationsTK.Lucon.PrintTrace()` for a convergence trace. See
  [`rotate`](@ref).
