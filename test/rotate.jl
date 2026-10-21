@testmodule NPLReference begin
    using DFTK
    using OrbitalRotationsTK

    # slow reference for the tests: ⟨φ_i|σ_F|φ_i⟩ of the NPL operators integrated on the
    # real-space grid, independent of the overlap densities used by OrbitalSubspace; the
    # Fourier coefficients of the potentials are shared with the production code
    function expectation_values(basis, φ)
        kpt = only(basis.kpoints)
        V = OrbitalRotationsTK.atom_local_potentials_real(basis)
        ρ = [abs2.(DFTK.ifft(basis, kpt, φi)) for φi in eachcol(φ)]
        [sum(ρi .* VF) * basis.dvol for ρi in ρ, VF in V]
    end
end


# The thresholds are about 10-100 times the run-to-run spread observed with an SCF tolerance
# of 1e-12 and a rotate tolerance of 1e-9; the reference values then depend only on the
# system, not on the random initial guess of the SCF or of rotate.

@testitem "Water ground state" tags=[:rotate] setup=[TestSystems] begin
    # a failure here, not in the tests below, points to a change in DFTK's SCF
    @test TestSystems.water.energies.total ≈ -17.5520470532783 atol=1e-10
end


@testitem "NPL loss at a fixed U" tags=[:rotate] setup=[TestSystems, NPLReference] begin
    using LinearAlgebra
    using Random

    basis = TestSystems.water.basis
    ψ = TestSystems.water_ψ
    U = OrbitalRotationsTK.random_unitary(Xoshiro(3), eltype(only(ψ)), size(only(ψ), 2))

    prep = OrbitalRotationsTK.prepare_gradient(NPL(), OrbitalSubspace(), basis, ψ)
    _, L = OrbitalRotationsTK.gradient(prep, [U], true)
    σ = NPLReference.expectation_values(basis, only(ψ) * U)
    @test L ≈ sum(abs2, σ) rtol=1e-13
end


@testitem "rotate water with NPL" tags=[:rotate] setup=[TestSystems, NPLReference] begin
    using LinearAlgebra
    using Random
    using OrbitalRotationsTK: Lucon

    basis = TestSystems.water.basis
    ψ = TestSystems.water_ψ
    res = rotate(basis, ψ, NPL(); rng=Xoshiro(1), tol=1e-9)
    U = res.optimizer.U
    φ = only(res.ψ)

    @test Lucon.converged(res.optimizer)
    @test U'U ≈ I atol=1e-12
    @test φ ≈ only(ψ) * U

    # reference loss of this system; guards against changes of the pseudopotential and
    # Fourier conventions in DFTK
    @test res.optimizer.loss ≈ 57.65136214439 rtol=1e-11

    σ = NPLReference.expectation_values(basis, φ)
    @test sum(abs2, σ) ≈ res.optimizer.loss rtol=1e-13

    # the mirror plane between the two H atoms maps the localized orbitals onto each other
    @test sort(σ[:, 2]) ≈ sort(σ[:, 3]) atol=1e-9

    # two O–H bonds: each H atom has one orbital with a large weight on it (about -0.63 Ha
    # versus at most -0.26 Ha for the others); the other two orbitals are lone pairs on O
    @test count(<(-0.5), σ[:, 2]) == count(<(-0.5), σ[:, 3]) == 1
    @test findfirst(<(-0.5), σ[:, 2]) != findfirst(<(-0.5), σ[:, 3])

    # any random start reaches the same optimum
    res2 = rotate(basis, ψ, NPL(); rng=Xoshiro(2), tol=1e-9)
    @test res2.optimizer.loss ≈ res.optimizer.loss rtol=1e-12
end
