using Documenter

makedocs(
    sitename = "ComputationalHeatTransfer",
    format = Documenter.HTML(),
    # modules = [ComputationalHeatTransfer],
    pages = [
        "Home" => "index.md",
        "Simulation" => "OHP simulation/OHP simulation.md",
        "Post-Processing" => "PostProcessing-oneresult/PostProcessing-oneresult.md",
        "OHP-DIY" => "OHP DIY/OHP DIY.md"
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "<github.com/liyuxuan48/ComputationalHeatTransfer.jl.git>",
    target = "build",
    deps = nothing,
    make = nothing,
    versions = nothing
)

