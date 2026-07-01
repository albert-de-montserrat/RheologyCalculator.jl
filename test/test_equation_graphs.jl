function equation_graph(c)
    return map(generate_equations(c)) do eq
        (;
            self = eq.self,
            parent = eq.parent,
            child = eq.child,
            fn = nameof(eq.fn),
            ind_input = eq.ind_input,
            el_number = eq.el_number,
        )
    end
end

@testset "equation graphs for composite layouts" begin
    v1 = LinearViscosity(1.0)
    v2 = LinearViscosity(2.0)
    v3 = LinearViscosity(3.0)
    v4 = LinearViscosity(4.0)
    v5 = LinearViscosity(5.0)
    e1 = IncompressibleElasticity(6.0)
    e2 = Elasticity(7.0, 8.0)
    dp = DruckerPrager(9.0, 10.0, 0.0)

    cases = (
        SeriesModel(v1, v2) => (
            (; self = 1, parent = 0, child = (), fn = :compute_strain_rate, ind_input = 1, el_number = (1, 2)),
        ),
        SeriesModel(v1, e2) => (
            (; self = 1, parent = 0, child = (), fn = :compute_strain_rate, ind_input = 1, el_number = (1, 1)),
            (; self = 2, parent = 0, child = (), fn = :compute_volumetric_strain_rate, ind_input = 2, el_number = (1, 1)),
        ),
        SeriesModel(v1, dp) => (
            (; self = 1, parent = 0, child = (2,), fn = :compute_strain_rate, ind_input = 1, el_number = (1, 1)),
            (; self = 2, parent = 1, child = (), fn = :compute_lambda, ind_input = 0, el_number = (1, 1)),
        ),
        SeriesModel(v1, e2, dp) => (
            (; self = 1, parent = 0, child = (2,), fn = :compute_strain_rate, ind_input = 1, el_number = (1, 1, 1)),
            (; self = 2, parent = 1, child = (), fn = :compute_lambda, ind_input = 0, el_number = (1, 1, 1)),
            (; self = 3, parent = 0, child = (), fn = :compute_volumetric_strain_rate, ind_input = 2, el_number = (1, 1, 1)),
        ),
        SeriesModel(ParallelModel(v1, v2)) => (
            (; self = 1, parent = 0, child = (2,), fn = :compute_strain_rate, ind_input = 1, el_number = ()),
            (; self = 2, parent = 1, child = (), fn = :compute_stress, ind_input = 0, el_number = (1, 2)),
        ),
        SeriesModel(v1, ParallelModel(v2, v3)) => (
            (; self = 1, parent = 0, child = (2,), fn = :compute_strain_rate, ind_input = 1, el_number = (1,)),
            (; self = 2, parent = 1, child = (), fn = :compute_stress, ind_input = 0, el_number = (2, 3)),
        ),
        SeriesModel(v1, ParallelModel(v2, SeriesModel(v3, e1))) => (
            (; self = 1, parent = 0, child = (2,), fn = :compute_strain_rate, ind_input = 1, el_number = (1,)),
            (; self = 2, parent = 1, child = (3,), fn = :compute_stress, ind_input = 0, el_number = (2,)),
            (; self = 3, parent = 2, child = (), fn = :compute_strain_rate, ind_input = 0, el_number = (3, 1)),
        ),
        SeriesModel(v1, ParallelModel(v2, e2)) => (
            (; self = 1, parent = 0, child = (2,), fn = :compute_strain_rate, ind_input = 1, el_number = (1,)),
            (; self = 2, parent = 1, child = (), fn = :compute_stress, ind_input = 0, el_number = (2, 1)),
            (; self = 3, parent = 0, child = (4,), fn = :compute_volumetric_strain_rate, ind_input = 2, el_number = (1,)),
            (; self = 4, parent = 3, child = (), fn = :compute_pressure, ind_input = 0, el_number = (2, 1)),
        ),
        SeriesModel(v1, ParallelModel(v2, dp)) => (
            (; self = 1, parent = 0, child = (2,), fn = :compute_strain_rate, ind_input = 1, el_number = (1,)),
            (; self = 2, parent = 1, child = (3, 4), fn = :compute_stress, ind_input = 0, el_number = (2, 1)),
            (; self = 3, parent = 2, child = (), fn = :compute_lambda_parallel, ind_input = 0, el_number = (2, 1)),
            (; self = 4, parent = 2, child = (), fn = :compute_plastic_strain_rate, ind_input = 0, el_number = (2, 1)),
        ),
        SeriesModel(v1, ParallelModel(v2, v3), ParallelModel(v4, SeriesModel(v5, e1))) => (
            (; self = 1, parent = 0, child = (2, 3), fn = :compute_strain_rate, ind_input = 1, el_number = (1,)),
            (; self = 2, parent = 1, child = (), fn = :compute_stress, ind_input = 0, el_number = (2, 3)),
            (; self = 3, parent = 1, child = (4,), fn = :compute_stress, ind_input = 0, el_number = (4,)),
            (; self = 4, parent = 3, child = (), fn = :compute_strain_rate, ind_input = 0, el_number = (5, 1)),
        ),
    )

    for (model, expected) in cases
        @test equation_graph(model) == expected
    end
end
