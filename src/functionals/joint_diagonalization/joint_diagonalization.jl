@doc raw"""
    JointDiagonalizationFunctional

Loss functionals of the form
```math
L(U) = \sum_F w_F \sum_i h\big(\langle \phi_i | \sigma_F | \phi_i \rangle\big),
\qquad \phi_i = \sum_j U_{ji} \, \psi_j,
```
with weights ``w_F``, Hermitian one-body operators ``\sigma_F`` and a scalar function
``h``. The Euclidean gradient is
```math
\Gamma_{pq} = \frac{\partial L}{\partial \overline{U_{pq}}}
= \sum_F w_F \, h'\big(\langle \phi_q | \sigma_F | \phi_q \rangle\big)
  \langle \psi_p | \sigma_F | \phi_q \rangle.
```

The loss is maximized, and ``h`` can be any scalar function. The family is named after the
case of a convex ``h`` (e.g. `Monomial(2)`): then the maximum is reached where all
``\sigma_F`` are as diagonal as possible in the rotated orbitals, i.e. where ``U`` jointly
diagonalizes them. For a non-convex ``h`` this interpretation no longer holds.

A subtype has the fields `h` (e.g. a [`Monomial`](@ref)) and `w` (one number for all ``F``,
or a vector with one weight per ``F``), and implements [`one_body_operators`](@ref) for
[`FourierSpace`](@ref) and [`RealSpace`](@ref). [`OrbitalSubspace`](@ref) is built from the
`FourierSpace` operators, so all three representations then work for it automatically.
"""
abstract type JointDiagonalizationFunctional end


@doc raw"""
    one_body_operators(functional, basis, ::RealSpace)
    one_body_operators(functional, basis, ::FourierSpace, Gs)

The one-body operators ``\sigma_F`` that define the `functional`: on the real-space grid,
or as an ``N_G \times N_F`` matrix of Fourier coefficients ``\sigma_F(\bm G)`` at the
vectors `Gs` (reduced coordinates).
"""
function one_body_operators end


maximize(::JointDiagonalizationFunctional) = true

# the expectation value ⟨ϕ_i|σ_F|ϕ_i⟩ is quadratic in U -> factor 2
max_taylor_degree(f::JointDiagonalizationFunctional) = 2 * taylor_degree(f.h)

# the weights w_F as a vector, from one number for all F or one number per F
function weights(f::JointDiagonalizationFunctional, NF::Integer)
    f.w isa Number && return fill(float(f.w), NF)
    length(f.w) == NF || throw(DimensionMismatch(
        "expected one weight per one-body operator ($NF), got $(length(f.w))"))
    return float.(collect(f.w))
end


@doc raw"""
    JointDiagOrbitalSubspaceCache

What [`gradient`](@ref) needs for a [`JointDiagonalizationFunctional`](@ref) in the
[`OrbitalSubspace`](@ref) representation, created by [`prepare_gradient`](@ref).

``\psi`` are the given orbitals and ``\phi`` the rotated ones. The orbitals are numbered
by a combined index ``(m, k)`` of orbital ``m`` and k-point ``k``, with the orbital running
fastest, over all k-points of the basis (spin included); at the Γ point there is a single
k-point.

- `h`, `w`: the scalar function ``h`` and the weights ``w_F``.
- `σ_ψψ`: ``\langle \psi_{mk} | \sigma_F | \psi_{nk'} \rangle`` as
  `σ_ψψ[(m,k), (n,k'), F]`.
- `σ_ψϕ`: scratch for ``\langle \psi_{pk} | \sigma_F | \phi_{qk'} \rangle``, same layout.
- `σ_ϕϕ_diag`: scratch for ``\langle \phi_q | \sigma_F | \phi_q \rangle`` as
  `σ_ϕϕ_diag[q, F]`.
- `Γ`: scratch for the gradient, one matrix per k-point, returned by `gradient`.
"""
struct JointDiagOrbitalSubspaceCache{H,TW,TS,TD,TG}
    h::H
    w::TW
    σ_ψψ::TS
    σ_ψϕ::TS
    σ_ϕϕ_diag::TD
    Γ::TG
end


function prepare_gradient(
    functional::JointDiagonalizationFunctional,
    representation::OrbitalSubspace,
    basis,
    ψ,
)
    # compute ⟨ψ_mk|e^{-iGr}|ψ_nk'⟩
    ρ, Gs = PsiTK.compute_overlap_densities(
        basis,
        ψ;
        callback=PsiTK.ShowProgress(desc="compute overlap densities"),
        representation.Ecut_ratio,
    )

    # compute Fourier representation of nuclear (pseudo-)potentials 
    V = one_body_operators(functional, basis, FourierSpace(), Gs)
    
    Nk, N, _, _, NG = size(ρ)
    NF = size(V, 2)

    # compute ⟨ψ_mk|σ^F|ψ_nk'⟩ 
    ρ_matrix = reshape(ρ, :, NG)
    # σ_ψψ[k,m,k',n,F] = Σ_G ρ[k,m,k',n,G] V[G,F]^*
    σ_ψψ = reshape(ρ_matrix * conj(V), Nk, N, Nk, N, NF)

    # switch index order k ⟷ m and k' ⟷ n, i.e. (1,2,3,4,5) ⟶ (2,1,4,3,5)
    σ_ψψ = permutedims(σ_ψψ, (2,1,4,3,5))
    σ_ψψ = reshape(σ_ψψ, N*Nk, N*Nk, NF)

    return JointDiagOrbitalSubspaceCache(
        functional.h,
        weights(functional, NF),
        σ_ψψ,
        similar(σ_ψψ),
        similar(σ_ψψ, real(eltype(σ_ψψ)), N, NF),
        [similar(σ_ψψ, N, N) for _ in 1:Nk],
    )
end


function gradient(prep::JointDiagOrbitalSubspaceCache, U, calc_loss)
    (; h, w, σ_ψψ, σ_ψϕ, σ_ϕϕ_diag, Γ) = prep
    rotate_kets!(σ_ψϕ, σ_ψψ, U)
    expectation_values!(σ_ϕϕ_diag, U, σ_ψϕ)
    euclidean_gradient!(Γ, h, w, σ_ϕϕ_diag, σ_ψϕ)
    return Γ, calc_loss ? loss(h, w, σ_ϕϕ_diag) : NaN
end


# the orbitals of k-point k in the combined index (orbital, k-point)
orbitals(N, k) = (k - 1) * N + 1:k * N


@doc raw"""
    rotate_kets!(σ_ψϕ, σ_ψψ, U)

Rotate the ket orbitals of the matrix elements in `σ_ψψ` with the unitaries `U`:
```math
\langle \psi_{pk} | \sigma_F | \phi_{qk'} \rangle
= \sum_n \langle \psi_{pk} | \sigma_F | \psi_{nk'} \rangle \, U^{(k')}_{nq}.
```
"""
function rotate_kets!(σ_ψϕ, σ_ψψ, U)
    N = size(first(U), 2)
    for F in axes(σ_ψψ, 3), (k′, U_k′) in enumerate(U)
        @views mul!(σ_ψϕ[:, orbitals(N, k′), F], σ_ψψ[:, orbitals(N, k′), F], U_k′)
    end
    return σ_ψϕ
end


@doc raw"""
    expectation_values!(σ_ϕϕ_diag, U, σ_ψϕ)

The expectation values of the operators in the rotated orbitals,
```math
\langle \phi_q | \sigma_F | \phi_q \rangle
= \sum_{k,k'} \sum_p \overline{U^{(k)}_{pq}} \,
  \langle \psi_{pk} | \sigma_F | \phi_{qk'} \rangle,
```
which are real, since ``\sigma_F`` is Hermitian. For several k-points, the sums still have
to be restricted to the k-points of one spin channel and normalized; so far only the Γ point
is supported.
"""
function expectation_values!(σ_ϕϕ_diag, U, σ_ψϕ)
    N = size(σ_ϕϕ_diag, 1)
    fill!(σ_ϕϕ_diag, 0)
    for F in axes(σ_ψϕ, 3), k′ in eachindex(U), (k, U_k) in enumerate(U), q in 1:N
        column = orbitals(N, k′)[q]
        @views σ_ϕϕ_diag[q, F] += real(dot(U_k[:, q], σ_ψϕ[orbitals(N, k), column, F]))
    end
    return σ_ϕϕ_diag
end


@doc raw"""
    euclidean_gradient!(Γ, h, w, σ_ϕϕ_diag, σ_ψϕ)

The Euclidean gradient of the loss, one matrix per k-point:
```math
\Gamma^{(k)}_{pq} = \frac{\partial L}{\partial \overline{U^{(k)}_{pq}}}
= \sum_F w_F \, h'\big(\langle \phi_q | \sigma_F | \phi_q \rangle\big)
  \sum_{k'} \langle \psi_{pk} | \sigma_F | \phi_{qk'} \rangle.
```
"""
function euclidean_gradient!(Γ, h, w, σ_ϕϕ_diag, σ_ψϕ)
    N = size(σ_ϕϕ_diag, 1)
    for (k, Γ_k) in enumerate(Γ)
        fill!(Γ_k, 0)
        for F in eachindex(w), k′ in eachindex(Γ), q in 1:N
            c = w[F] * derivative(h, σ_ϕϕ_diag[q, F])
            column = orbitals(N, k′)[q]
            @views Γ_k[:, q] .+= c .* σ_ψϕ[orbitals(N, k), column, F]
        end
    end
    return Γ
end


@doc raw"""
    loss(h, w, σ_ϕϕ_diag)

The loss
```math
L = \sum_F w_F \sum_q h\big(\langle \phi_q | \sigma_F | \phi_q \rangle\big).
```
"""
loss(h, w, σ_ϕϕ_diag) =
    sum(w[F] * h(σ_ϕϕ_diag[q, F]) for q in axes(σ_ϕϕ_diag, 1), F in eachindex(w))
