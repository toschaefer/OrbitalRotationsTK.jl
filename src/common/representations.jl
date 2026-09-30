@doc raw"""
    OrbitalSubspace(; Ecut_ratio=4.0)

Evaluate loss functionals via the matrix elements
``\langle \psi_m | \sigma_F | \psi_n \rangle`` of their operators ``\sigma_F`` in the
subspace spanned by the orbitals ``\psi``, computed from the overlap densities of the
orbitals.

`Ecut_ratio` is the ratio of the plane-wave cutoff of the overlap densities to the orbital
cutoff. Values up to `supersampling^2` of the basis (4 for DFTK's default) are allowed,
since the FFT grid holds products of orbitals exactly up to there; `Ecut_ratio=4` therefore
yields the exact overlap densities.

See [Functionals and representations](@ref) for the available combinations and their cost.
"""
@kwdef struct OrbitalSubspace
    Ecut_ratio::Float64 = 4.0
end


"""
    FourierSpace(; Ecut_ratio=4.0)

Evaluate loss functionals using Fourier representation of integrals. Not yet implemented.

`Ecut_ratio` has the same meaning as in [`OrbitalSubspace`](@ref).

See [Functionals and representations](@ref) for the available combinations and their cost.
"""
@kwdef struct FourierSpace
    Ecut_ratio::Float64 = 4.0
end


"""
    RealSpace()

Evaluate loss functionals using real space representation of integrals. Not yet
implemented.

See [Functionals and representations](@ref) for the available combinations and their cost.
"""
struct RealSpace end


"""
    prepare_gradient(functional, representation, basis, ψ)

Precompute and preallocate what the `functional` needs for the given `representation`, for
the orbitals `ψ` of `basis` (one matrix per k-point). Runs once, before the optimization;
the result is passed to [`gradient`](@ref) in every step.
"""
function prepare_gradient end


@doc raw"""
    gradient(prep, U, calc_loss) -> (Γ, L)

Loss `L` and Euclidean gradient `Γ` at the unitaries `U`, using `prep` from
[`prepare_gradient`](@ref). `U` and `Γ` hold one matrix per k-point ``k``, in the order of
`basis.kpoints` (spin included), with
```math
\Gamma^{(k)}_{pq} = \frac{\partial L}{\partial \overline{U^{(k)}_{pq}}}.
```
`L` is only computed if `calc_loss` is true, and is `NaN` otherwise. `Γ` may be a buffer
that the next call overwrites. Called by the optimizer in every step.
"""
function gradient end


"""
    maximize(functional)

Whether the loss of `functional` is maximized (`true`) or minimized (`false`). Each family
of functionals implements it.
"""
function maximize end


"""
    max_taylor_degree(functional)

The polynomial degree of the loss of `functional` in the entries of the unitary, or an
effective degree for a non-polynomial loss (see [`taylor_degree`](@ref)). The line search
of the optimizer uses it to choose its step sizes. Each family of functionals implements
it.
"""
function max_taylor_degree end
