using RheologyCalculator, Test
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic
import RheologyCalculator: compute_residual

include("../rheologies/RheologyDefinitions.jl")
include("../examples/tensor_helpers.jl")

function runtests()
    files = readdir(@__DIR__)
    test_files = filter(startswith("test_"), files)

    allocations_only = "--allocations-only" in ARGS
    if allocations_only
        test_files = filter(==("test_allocations.jl"), test_files)
    elseif Base.JLOptions().code_coverage != 0
        filter!(!=("test_allocations.jl"), test_files)
    end

    for f in test_files
        !isdir(f) && include(f)
    end
    return
end

runtests()
