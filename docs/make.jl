push!(LOAD_PATH, "../src/")

using Documenter
using PBCCompiler

DocMeta.setdocmeta!(PBCCompiler, :DocTestSetup, :(using QuantumClifford, PBCCompiler); recursive=true)

makedocs(
    doctest = false,
    clean = true,
    warnonly = :missing_docs,
    sitename = "PBCCompiler.jl",
    format = Documenter.HTML(),
    modules = [PBCCompiler, PBCCompiler.Verification, PBCCompiler.GadgetShapes, PBCCompiler.ResourceEstimation],
    authors = "Stefan Krastanov",
    pages = [
        "PBCCompiler.jl" => "index.md",
        "Verification" => "verification.md",
        "Resource estimation" => "resource_estimation.md",
        "API" => "API.md",
    ]
)

deploydocs(
    repo = "github.com/QuantumSavory/PBCCompiler.jl.git",
    push_preview = true,
)
