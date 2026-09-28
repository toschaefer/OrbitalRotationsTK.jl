"""
Result of [`rotate`](@ref): the rotated orbitals `ψ` (one matrix per k-point), the result
of the unitary optimization `optimizer` (with the optimal unitary, loss and convergence
information, see `Lucon.optimize`), and the `functional` and `representation` used.
"""
struct RotationResult{Tψ,TO,TF,TR}
    ψ::Tψ
    optimizer::TO
    functional::TF
    representation::TR
end


"""
    rotate(basis, ψ, functional; representation=OrbitalSubspace(), U0=nothing,
           tol=1e-6, maxiter=1000, callback=nothing)

Find the unitary `U` that optimizes the `functional` for the rotated orbitals
ϕ_i = Σ_j U_ji ψ_j, i.e. `ψk * U` for the orbitals `ψk` of a k-point, and return a
[`RotationResult`](@ref). Currently Γ-point only.

- `basis`: the `PlaneWaveBasis` of the orbitals, e.g. `scfres.basis`.
- `ψ`: orthonormal orbitals of `basis`, one matrix per k-point with one orbital per column,
  as in `scfres.ψ`.
- `representation`: how the loss and its gradient are evaluated.
- `U0`: initial unitaries, one per k-point; the identity if `nothing`.
- `tol`: convergence threshold for the largest element |G_pq| of the Riemannian gradient
  G = ΓU† − UΓ† of the loss, Γ_pq = ∂L/∂conj(U_pq), in units of the loss. G_pq is
  proportional to the derivative of the loss with respect to the angle of a rotation of
  the orbital pair (p,q), so the criterion does not grow with the number of orbitals. The
  units of the loss are given in the docstring of each functional.
- `maxiter`: maximum number of rotations of `U`.
- `callback`: called once per iteration with the state of the optimizer; returning `true`
  stops the optimization.

The optimization is done by `Lucon.optimize`; see its documentation for the algorithm, the
state passed to the callback, and ready-made callbacks such as a convergence trace.
"""
function rotate(
    basis,
    ψ::AbstractVector{<:AbstractMatrix},
    functional;
    representation=OrbitalSubspace(),
    U0::Union{Nothing,AbstractVector{<:AbstractMatrix}}=nothing,
    tol=1e-6,
    maxiter=1000,
    callback=nothing,
)
    # consistency checks
    length(basis.kpoints) == 1 && iszero(basis.kpoints[1].coordinate) ||
        throw(ArgumentError("currently Γ-only!"))
    length(ψ) == length(basis.kpoints) || throw(DimensionMismatch(
        "expected orbitals for $(length(basis.kpoints)) k-point(s), got $(length(ψ))"))
    ψk = only(ψ)
    N = size(ψk, 2)
    U = isnothing(U0) ? Matrix{eltype(ψk)}(I, N, N) : only(U0)
    size(U) == (N, N) || throw(DimensionMismatch("U0 must be $N×$N, got $(size(U))"))
    U'U ≈ I || throw(ArgumentError("U0 must be unitary"))

    # prepare and optimize
    prep = prepare_gradient(functional, representation, basis, ψ)
    res = Lucon.optimize(
        (U, calc_loss) -> gradient(prep, U, calc_loss),
        U;
        max_taylor_degree=max_taylor_degree(functional),
        maximize=maximize(functional),
        max_gradient_tolerance=tol,
        max_iter=maxiter,
        callback,
    )
    return RotationResult([ψk * res.U], res, functional, representation)
end
