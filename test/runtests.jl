using RheologyCalculator, Test
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic
import RheologyCalculator: compute_residual

include("../examples/RheologyDefinitions.jl")

function runtests()
    files = readdir(@__DIR__)
    test_files = filter(startswith("test_"), files)

    for f in test_files
        !isdir(f) && include(f)
    end
    return
end

runtests()