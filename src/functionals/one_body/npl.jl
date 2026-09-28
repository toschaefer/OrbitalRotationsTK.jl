"""
    NPL(; h=Monomial(2), w=1.0)

The Nuclear Potential Localization (NPL) functional

    L(U) = Σ_F w_F Σ_i h(⟨ϕ_i|σ_F|ϕ_i⟩),    ϕ_i = Σ_j U_ji ψ_j

The one-body operators σ_F are the local (attractive) pseudopotentials of the atoms F,
hence ⟨ϕ_i|σ_F|ϕ_i⟩ < 0 and `h` must be defined for negative arguments. The expectation
values are energies in Hartree, so with `Monomial(p)` the loss has units of Hartreeᵖ
(times the units of `w`), which sets the scale of `tol` in [`rotate`](@ref).

The weights `w` are either one number for all atoms or a vector with one weight per atom, in
the order of `basis.model.atoms`. A weight of zero excludes the atom from the functional.
"""
@kwdef struct NPL{H,TW} <: OneBodyFunctional
    h::H = Monomial(2)
    w::TW = 1.0
end


function one_body_operators(::NPL, basis, ::FourierSpace, Gs)
    # - N_G × N_F matrix of the atomic pseudopotentials σ_F(G) at the G vectors Gs:
    #   reduce(hcat, atom_local_potentials_fourier(basis, Gs))
end



# Local pseudopotential of each atom in Fourier space, evaluated at the G vectors `Gs`
# (fractional coordinates). Returns one coefficient vector per atom, ordered like `Gs`.
#
# Uses the factorization v_A(G) = v̂_A(|G|) * cis2pi(-G⋅R_A) / sqrt(Ω) into a radial form
# factor and a structure factor. Evaluating v̂ is a radial quadrature, so it is done once per
# distinct |G| and element instead of once per (G, atom) pair.
function atom_local_potentials_fourier(basis, Gs)
    model = basis.model
    Ω = model.unit_cell_volume
    T = eltype(model.recip_lattice)

    # one radial quadrature per distinct |G| and per element group
    ip   = Vector{Int}(undef, length(Gs))   # index into Gs -> index into ps
    ps   = T[]
    seen = Dict{T,Int}()
    for (i, G) in enumerate(Gs)
        p = norm(model.recip_lattice * G)
        ip[i] = get!(seen, p) do
            push!(ps, p); length(ps)
        end
    end
    form_factors = [DFTK.local_potential_fourier(model.atoms[first(group)], ps)
                    for group in model.atom_groups]

    # the form factors follow model.atom_groups, so invert atom -> group
    igroup_of_atom = zeros(Int, length(model.atoms))
    for (igroup, group) in enumerate(model.atom_groups), iatom in group
        igroup_of_atom[iatom] = igroup
    end

    map(eachindex(model.atoms)) do iatom
        R  = model.positions[iatom]
        ff = form_factors[igroup_of_atom[iatom]]
        [ff[ip[i]] * DFTK.cis2pi(-G⋅R) / sqrt(Ω) for (i, G) in enumerate(Gs)]
    end
end


# Real-space local pseudopotential of each atom, one `basis.fft_size` array per atom.
# Needs the full cubic FFT grid, since `irfft` consumes every coefficient.
atom_local_potentials_real(basis) =
    [DFTK.irfft(basis, reshape(V, basis.fft_size))
     for V in atom_local_potentials_fourier(basis, vec(G_vectors(basis)))]

