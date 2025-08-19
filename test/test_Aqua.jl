using Aqua, Test, RheologyCalculator

@testset "Project extras" begin
    @test Aqua.test_project_extras(RheologyCalculator).value
end

@testset "Undefined exports" begin
    @test Aqua.test_undefined_exports(RheologyCalculator).value
end

@testset "Compats" begin
    @test !Aqua.test_deps_compat(
        RheologyCalculator;
        check_julia = true,
        check_extras = false,
    ).anynonpass
    # @test Aqua.test_stale_deps(RheologyCalculator).value
end

@testset "Persistent tasks" begin
    Aqua.test_persistent_tasks(RheologyCalculator)
end

@testset "Ambiguities" begin
    @test Aqua.test_ambiguities(RheologyCalculator).value
end

@testset "Piracy" begin
    @test Aqua.test_piracies(RheologyCalculator).value
end
