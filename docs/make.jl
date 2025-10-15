using Documenter, OceanLight

makedocs(
    sitename="OceanLight.jl", format = Documenter.HTML(prettyurls=true),
    pages = [
        "Home" => "index.md",
        "Quick Start" => [
        "Photons at the center" => "QuickStart/Center.md",
        "Implement MPI" => "QuickStart/MPI.md"
        ],
        "Module" => [
        "Parameters" => "module/Parameters.md",
        "Refraction" => "module/refraction.md",
        "Monte Carlo Simulation" => "module/MonteCarlo.md"
        ],
        "Simulation" => [
        "Model Setup" => "Simulation/ModelSetup.md",
        "Air-Water Interaction" => "Simulation/AirWaveInteract.md",
        "Light within Water" => "Simulation/WithinWater.md",
        "Exporting data" => "Simulation/Exporting.md"
        ],
        "Contributor's guide" => "contribute.md",
        "Reference" => "reference.md"
    ])

deploydocs(;
    repo="github.com/haoboatlab/OceanLight.jl",
)
