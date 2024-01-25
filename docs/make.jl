using Documenter, LightMC

makedocs(
    sitename="LightMC.jl", format = Documenter.HTML(prettyurls = false),
    pages = [
        "Home" => "index.md",
        "Quick Start" => [
        "Photons at the center" => "QuickStart/Center.md"
        "Everywhere on the water surface" => "QuickStart/Everywhere.md"
        ],
        "module" => [
        "Parameters" => "module/Parameters.md",
        "Refraction" => "module/refraction.md",
        "Monte Carlo Simulation" => "module/MonteCarlo.md"
        ],
        "Reference" => "reference.md"
    ])

deploydocs(;
    repo="github.com/haoboatlab/LightMC.jl",
)