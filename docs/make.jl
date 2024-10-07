using Documenter, OceanLight

makedocs(
    sitename="OceanLight.jl", format = Documenter.HTML(prettyurls = false),
    pages = [
        "Home" => "index.md",
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
        "Reference" => "reference.md"
    ])

deploydocs(;
    repo="github.com/haoboatlab/OceanLight.jl",
)