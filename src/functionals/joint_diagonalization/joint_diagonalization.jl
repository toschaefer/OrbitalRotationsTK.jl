"""
    JointDiagonalizationFunctional

Loss functionals of the form

    L(U) = Σ_F w_F Σ_i h(⟨ϕ_i|σ_F|ϕ_i⟩),    ϕ_i = Σ_j U_ji ψ_j

with weights w_F, Hermitian one-body operators σ_F and a scalar function h.
Euclidean gradient:

    Γ_pq = ∂L/∂conj(U_pq) = Σ_F w_F h'(⟨ϕ_q|σ_F|ϕ_q⟩) ⟨ψ_p|σ_F|ϕ_q⟩

The loss is maximized, and `h` can be any scalar function. The family is named after the
case of a convex `h` (e.g. `Monomial(2)`): then the maximum is reached where all σ_F are as
diagonal as possible in the rotated orbitals, i.e. where U jointly diagonalizes them. For a
non-convex `h` this interpretation no longer holds.

A subtype has the fields `h` (e.g. a `Monomial`) and `w` (one number for all F, or a vector
with one weight per F), and implements `one_body_operators` for `FourierSpace` and
`RealSpace`. `OrbitalSubspace` is built from the `FourierSpace` operators, so all three
representations then work for it automatically.
"""
abstract type JointDiagonalizationFunctional end


"""
    one_body_operators(functional, basis, ::RealSpace)
    one_body_operators(functional, basis, ::FourierSpace, Gs)

The one-body operators σ_F that define the `functional`: on the real-space grid, or as an
N_G × NF matrix of Fourier coefficients σ_F(G) at the G vectors `Gs`.
"""
function one_body_operators end


maximize(::JointDiagonalizationFunctional) = true

# the expectation value ⟨ϕ_i|σ_F|ϕ_i⟩ is quadratic in U -> factor 2
max_taylor_degree(f::JointDiagonalizationFunctional) = 2 * taylor_degree(f.h)

function weights(f::JointDiagonalizationFunctional, NF::Integer)
    f.w isa Number && return fill(float(f.w), NF)
    length(f.w) == NF || throw(DimensionMismatch(
        "expected one weight per one-body operator ($NF), got $(length(f.w))"))
    return float.(collect(f.w))
end


# ψ: the given orbitals, ϕ: the rotated ones. Combined index (orbital, k-point) with the
# orbital fastest, over all k-points of the basis (spin included); at Γ a single k-point
struct JointDiagOrbitalSubspaceCache{H,TW,TS,TD,TG}
    h::H
    w::TW             # weights w_F
    σ_ψψ::TS          # ⟨ψ_mk|σ_F|ψ_nk'⟩
    σ_ψϕ::TS          # scratch: ⟨ψ_pk|σ_F|ϕ_qk'⟩
    σ_ϕϕ_diag::TD     # scratch: ⟨ϕ_q|σ_F|ϕ_q⟩
    Γ::TG             # scratch: Γ[k] = ∂L/∂conj(U[k]), returned by gradient
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


# ⟨ψ_pk|σ_F|ϕ_qk'⟩ = Σ_n ⟨ψ_pk|σ_F|ψ_nk'⟩ U[k'][n,q]
function rotate_kets!(σ_ψϕ, σ_ψψ, U)
    N = size(first(U), 2)
    for F in axes(σ_ψψ, 3), (k′, U_k′) in enumerate(U)
        @views mul!(σ_ψϕ[:, orbitals(N, k′), F], σ_ψψ[:, orbitals(N, k′), F], U_k′)
    end
    return σ_ψϕ
end


# ⟨ϕ_q|σ_F|ϕ_q⟩ = Σ_{k,k'} Σ_p conj(U[k][p,q]) ⟨ψ_pk|σ_F|ϕ_qk'⟩, summed over the k-points
# of the spin channel of q (normalization for k-points to be fixed); real, since σ_F is
# Hermitian
function expectation_values!(σ_ϕϕ_diag, U, σ_ψϕ)
    N = size(σ_ϕϕ_diag, 1)
    fill!(σ_ϕϕ_diag, 0)
    for F in axes(σ_ψϕ, 3), k′ in eachindex(U), (k, U_k) in enumerate(U), q in 1:N
        column = orbitals(N, k′)[q]
        @views σ_ϕϕ_diag[q, F] += real(dot(U_k[:, q], σ_ψϕ[orbitals(N, k), column, F]))
    end
    return σ_ϕϕ_diag
end


# Γ[k][p,q] = ∂L/∂conj(U[k][p,q]) = Σ_F w_F h′(⟨ϕ_q|σ_F|ϕ_q⟩) Σ_k' ⟨ψ_pk|σ_F|ϕ_qk'⟩
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


# L = Σ_F w_F Σ_q h(⟨ϕ_q|σ_F|ϕ_q⟩)
loss(h, w, σ_ϕϕ_diag) =
    sum(w[F] * h(σ_ϕϕ_diag[q, F]) for q in axes(σ_ϕϕ_diag, 1), F in eachindex(w))
