using Documenter
using OrbitalRotationsTK

makedocs(;
    sitename="OrbitalRotationsTK.jl",
    modules=[OrbitalRotationsTK],
    repo=Remotes.GitHub("toschaefer", "OrbitalRotationsTK.jl"),
    format=Documenter.HTML(; edit_link="main"),
    checkdocs=:exports,
    pages=[
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Functionals and representations" => "functionals_and_representations.md",
        "Code reference" => "code_reference.md",
    ],
)

deploydocs(; repo="github.com/toschaefer/OrbitalRotationsTK.jl", devbranch="main")
