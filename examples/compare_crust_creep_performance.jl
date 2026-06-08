using GLMakie
using Printf
using Statistics

GLMakie.activate!(; visible = false)

module CrustRheologyCalculatorExample
include("Elastic_Diffusion_Dislocation_Crust.jl")
end

module CrustNewtonExample
include("Elastic_Diffusion_Dislocation_Crust_Newton.jl")
end

const RC = CrustRheologyCalculatorExample
const NR = CrustNewtonExample

function run_rheologycalculator(setup; ntime, dt)
    return RC.simulate_crustal_creep(
        setup.c,
        setup.diffusion,
        setup.dislocation,
        setup.vars,
        setup.x,
        setup.xnorm,
        setup.env;
        ntime = ntime,
        dt = dt,
    )
end

function run_newton(setup; ntime, dt)
    return NR.simulate_crustal_creep_newton(
        setup.diffusion,
        setup.dislocation,
        setup.elastic,
        setup.vars,
        setup.env;
        ntime = ntime,
        dt = dt,
    )
end

function measure(label, f; samples = 20, warmups = 3)
    for _ in 1:warmups
        f()
    end

    times = zeros(samples)
    bytes = zeros(Int, samples)
    result = nothing

    GC.gc()
    for i in 1:samples
        timed = @timed f()
        result = timed.value
        times[i] = timed.time
        bytes[i] = timed.bytes
    end

    return (; label, times, bytes, result)
end

function print_result(r)
    @printf(
        "%-28s %10.4f %10.4f %10.4f %14.2f\n",
        r.label,
        1.0e3 * minimum(r.times),
        1.0e3 * median(r.times),
        1.0e3 * maximum(r.times),
        median(r.bytes) / 1024^2,
    )
end

function compare_results(rc_history, nr_history)
    rc_τ = rc_history.τII
    nr_τ = nr_history.τII

    max_abs_τ = maximum(abs.(rc_τ .- nr_τ))
    max_rel_τ = maximum(abs.(rc_τ .- nr_τ) ./ max.(abs.(rc_τ), eps(Float64)))

    return (; max_abs_τ, max_rel_τ)
end

function main(; ntime = 501, dt = 250 * RC.SecYear, samples = 20, warmups = 3)
    rc_setup = RC.setup_crustal_creep()
    nr_setup = NR.setup_crustal_creep_newton()

    println("Comparing crustal elastic + diffusion + dislocation examples")
    println("ntime = $ntime, dt = $(dt / RC.SecYear) yr, samples = $samples, warmups = $warmups")
    println()

    rc = measure(
        "RheologyCalculator.jl",
        () -> run_rheologycalculator(rc_setup; ntime = ntime, dt = dt);
        samples = samples,
        warmups = warmups,
    )

    nr = measure(
        "manual Newton-Raphson",
        () -> run_newton(nr_setup; ntime = ntime, dt = dt);
        samples = samples,
        warmups = warmups,
    )

    println("timing excludes package/script load time and figure creation")
    @printf("%-28s %10s %10s %10s %14s\n", "case", "min ms", "median ms", "max ms", "median MiB")
    print_result(rc)
    print_result(nr)

    speedup = median(rc.times) / median(nr.times)
    consistency = compare_results(rc.result, nr.result)

    println()
    @printf("Newton median speedup: %.2fx\n", speedup)
    @printf("max |ΔτII|: %.6e Pa (%.6e relative)\n", consistency.max_abs_τ, consistency.max_rel_τ)

    return (; rheologycalculator = rc, newton = nr, speedup, consistency)
end

out = main();
