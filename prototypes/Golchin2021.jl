using GLMakie

yield(p, q, A, B, C, β) = (p-C)^2/A^2 + (q'-β*p)^2/B^2 - 1

Af(p, pc, pt, γ) = (pc - pt)/(2*π) *(2*atan(γ*(pc+pt-2p)/(2*pc))+π)

Bf(p, pc, pt, M, C, α) = M*C*exp(α*(p - C)/(pc - pt))

Cf(pc, pt, γ) = (pc - pt)/π * atan(γ/2) + (pc + pt)/2 

function main()

    pc  = 2e8
    p   = LinRange(-0.7e8, 2.1e8, 100)
    q   = LinRange(-0.0e8, 1.8e8, 100)

    pt  = -1e7
    β   = 0.

    γ = 0.5
    α = 0.5

    ###################################

    M  = 1.0 
    @show ϕ  = asind(3*M/(6+M))
    @show 6*sind(ϕ) / (3 - sind(ϕ))

    τcs = M*p

    C = Cf.(pc, pt, γ)
    B = Bf.(p, pc, pt, M, C, α)
    A = Af.(p, pc, pt, γ)
    f = yield.(p, q', A, B, C, β)
    @show extrema(f)

    ###################################

    M  = 0.6 
    @show ϕ  = asind(3*M/(6+M))
    @show 6*sind(ϕ) / (3 - sind(ϕ))

    C = Cf.(pc, pt, γ)
    B = Bf.(p, pc, pt, M, C, α)
    A = Af.(p, pc, pt, γ)
    g = yield.(p, q', A, B, C, β)
    @show extrema(g)

    ###################################

    fig = Figure()
    ax = Axis(fig[1,1], aspect=DataAspect())
    lines!(ax, p./pc, τcs./pc)
    contour!(ax, p./pc, q./pc, f, levels=[0])
    contour!(ax, p./pc, q./pc, g, levels=[0], linestyle=:dash)
    # ylims!(ax, 0, 0.8)
    display(fig)
end


main()