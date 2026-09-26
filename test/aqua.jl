@testitem "Aqua.jl quality assurance" tags=[:aqua] begin
    using Aqua
    Aqua.test_all(OrbitalRotationsTK)
end
