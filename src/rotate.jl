"""
Result of [`rotate`](@ref): the rotated orbitals `ψk`, the result of the unitary
optimization `optimizer` (with the optimal unitary, loss and convergence information, see
`Lucon.optimize`), and the `functional` and `representation` used.
"""
struct RotationResult{Tψk,TO,TF,TR}
    ψk::Tψk
    optimizer::TO
    functional::TF
    representation::TR
end


"""
    rotate(basis, ψk, functional; representation=OrbitalSubspace(), U0=nothing,
           tol=1e-6, maxiter=1000, callback=nothing)

Find the unitary `U` that optimizes the `functional` for the rotated orbitals
ϕ_i = Σ_j U_ji ψ_j, i.e. `ψk * U`, and return a [`RotationResult`](@ref).
Currently Γ-point only.

- `basis`: the `PlaneWaveBasis` of the orbitals, e.g. `scfres.basis`.
- `ψk`: orthonormal orbitals of `basis` at the Γ point, one per column.
- `representation`: how the loss and its gradient are evaluated.
- `U0`: initial unitary; the identity if `nothing`.
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
    ψk,
    functional;
    representation=OrbitalSubspace(),
    U0=nothing,
    tol=1e-6,
    maxiter=1000,
    callback=nothing,
)
    length(basis.kpoints) == 1 && iszero(basis.kpoints[1].coordinate) ||
        throw(ArgumentError("currently Γ-only!"))
    N = size(ψk, 2)
    U = @something U0 Matrix{eltype(ψk)}(I, N, N)
    size(U) == (N, N) || throw(DimensionMismatch("U0 must be $N×$N, got $(size(U))"))
    U'U ≈ I || throw(ArgumentError("U0 must be unitary"))

    prep = prepare_gradient(functional, representation, basis, ψk)

    res = Lucon.optimize(
        (U, calc_loss) -> gradient(prep, U, calc_loss),
        U;
        max_taylor_degree=max_taylor_degree(functional),
        maximize=maximize(functional),
        max_gradient_tolerance=tol,
        max_iter=maxiter,
        callback,
    )

    return RotationResult(ψk * res.U, res, functional, representation)
end
