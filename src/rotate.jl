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
    rotate(basis, ψk, functional; representation=OrbitalSubspace(), U0=nothing, kwargs...)

Find the unitary `U` that optimizes the `functional` for the rotated orbitals
ϕ_i = Σ_j U_ji ψ_j, i.e. `ψk * U`, and return a [`RotationResult`](@ref).
Currently Γ-point only.

- `basis`: the `PlaneWaveBasis` of the orbitals, e.g. `scfres.basis`.
- `ψk`: orthonormal orbitals of `basis` at the Γ point, one per column.
- `representation`: how the loss and its gradient are evaluated.
- `U0`: initial unitary; the identity if `nothing`.

All other keyword arguments, e.g. for convergence control or a callback, are passed to
`Lucon.optimize`; see its docstring.
"""
function rotate(
    basis,
    ψk,
    functional;
    representation=OrbitalSubspace(),
    U0=nothing,
    kwargs...,
)
    length(basis.kpoints) == 1 && iszero(basis.kpoints[1].coordinate) ||
        throw(ArgumentError("currently Γ-only!"))
    N = size(ψk, 2)
    U = @something U0 Matrix{eltype(ψk)}(I, N, N)
    size(U) == (N, N) || throw(DimensionMismatch("U0 must be $N×$N, got $(size(U))"))
    U'U ≈ I || throw(ArgumentError("U0 must be unitary"))

    prep = prepare_gradient(functional, representation, basis, ψk)

    # kwargs first: the functional fixes the Taylor degree and the direction
    res = Lucon.optimize(
        (U, calc_loss) -> gradient(prep, U, calc_loss),
        U;
        kwargs...,
        max_taylor_degree=max_taylor_degree(functional),
        maximize=maximize(functional),
    )

    return RotationResult(ψk * res.U, res, functional, representation)
end
