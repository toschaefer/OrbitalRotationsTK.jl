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
N_G × N_F matrix of Fourier coefficients σ_F(G) at the G vectors `Gs`.
"""
function one_body_operators end


maximize(::JointDiagonalizationFunctional) = true

# the expectation value ⟨ϕ_i|σ_F|ϕ_i⟩ is quadratic in U -> factor 2
max_taylor_degree(f::JointDiagonalizationFunctional) = 2 * taylor_degree(f.h)

function weights(f::JointDiagonalizationFunctional, n_F::Integer)
    f.w isa Number && return fill(float(f.w), n_F)
    length(f.w) == n_F || throw(DimensionMismatch(
        "expected one weight per one-body operator ($n_F), got $(length(f.w))"))
    return float.(collect(f.w))
end


struct JointDiagOrbitalSubspaceCache{H,TW,TS,TB}
    h::H
    w::TW          # weights w_F
    σ::TS          # N_F x N x N (⟨ψ_i|σ_F|ψ_j⟩)
    σ_rotated::TS  # scratch
    buffer_Fj::TB  # scratch
end


function prepare_gradient(
    functional::JointDiagonalizationFunctional,
    representation::OrbitalSubspace,
    basis,
    ψk
)
    # - ρ_ij(G), Gs: overlap densities of ψk up to representation.Ecut_ratio (from PsiTK,
    #   which needs a (basis, ψk) entry point)
    # - V = one_body_operators(functional, basis, FourierSpace(), Gs)      # N_G × N_F
    # - σ[F,i,j] = Σ_G conj(V[G,F]) ρ_ij(G); exact with DFTK's FFT normalization, no dvol
    # - w = weights(functional, size(V, 2))
    # - return JointDiagOrbitalSubspaceCache(functional.h, w, σ, similar(σ),
    #       similar(σ, N_F, N))
    
    ρmnG, Gs = PsiTK.compute_overlap_densities(
        basis,
        ψk; 
        callback = PsiTK.ShowProgress(desc="compute overlap densities"),
        Ecut_ratio = 4.0
    )
end


function gradient(prep::JointDiagOrbitalSubspaceCache, U, calc_loss)
    # - rotate the ket index into prep.σ_rotated: B[F,p,q] = Σ_j σ[F,p,j] U[j,q] (one gemm)
    # - diagonal into prep.buffer_Fj:  d[F,q] = Σ_p conj(U[p,q]) B[F,p,q] = ⟨ϕ_q|σ_F|ϕ_q⟩
    #   (real, since σ_F is Hermitian)
    # - Γ[p,q] = Σ_F w_F h′(d[F,q]) B[F,p,q]
    # - L = Σ_F w_F Σ_q h(d[F,q]) if calc_loss, else missing
    # - return (Γ, L)
end
