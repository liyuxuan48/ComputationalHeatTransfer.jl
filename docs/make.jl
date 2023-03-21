using Documenter

makedocs(
    sitename = "ComputationalHeatTransfer",
    format = Documenter.HTML(),
    # modules = [ComputationalHeatTransfer],
    pages = [
        "Home" => "index.md",
        "Simulation" => "manual/OHP simulation.md",
        "Post-Processing" => "manual/PostProcessing-oneresult.md",
        "OHP-DIY" => "manual/OHP DIY.md"
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

