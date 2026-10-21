@testmodule TestSystems begin
    using DFTK
    using LinearAlgebra
    using PseudoPotentialData

    # water in a cubic box at the Γ point, the system of the tutorial
    function water_scfres()
        pseudopotentials = PseudoFamily("dojo.nc.sr.lda.v0_4_1.standard.upf")
        O = ElementPsp(:O, pseudopotentials)
        H = ElementPsp(:H, pseudopotentials)

        a = 12.0
        lattice = a * I(3)
        positions_cart = [[0.0, 0.0, 0.0], [1.43, 1.11, 0.0], [-1.43, 1.11, 0.0]]
        positions = [r / a .+ 0.5 for r in positions_cart]

        model = model_DFT(lattice, [O, H, H], positions; functionals=LDA())
        basis = PlaneWaveBasis(model; Ecut=20, kgrid=(1, 1, 1))
        # tight, since the error of the orbitals enters the loss linearly
        self_consistent_field(basis; tol=1e-12, callback=identity)
    end

    const water = water_scfres()
    const water_ψ =
        DFTK.select_occupied_orbitals(water.basis, water.ψ, water.occupation).ψ
end
