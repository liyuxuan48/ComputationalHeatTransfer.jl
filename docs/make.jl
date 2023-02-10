using Documenter

makedocs(
    prettyurls = true,
    sitename = "ComputationalHeatTransfer",
    format = Documenter.HTML(),
    # modules = [ComputationalHeatTransfer],
    pages = [
        "Home" => "index.md",
        "Basics" => ["manual/OHP simulation.md"
                     ]
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "<github.com/liyuxuan48/OscillatingHeatPipe.git>",
    target = "build",
    deps = nothing,
    make = nothing,
    versions = nothing
)

