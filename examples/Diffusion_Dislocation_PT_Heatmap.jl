using RheologyCalculator

using GLMakie
using LaTeXStrings

include("../rheologies/RheologyDefinitions.jl")

function creep_prefactors(diffusion, dislocation, T, P, f, d)
    Cdiff = diffusion.A * f^diffusion.r * d^(-diffusion.p) * exp(-(diffusion.E + P * diffusion.V) / (diffusion.R * T))
    Cdisl = dislocation.A * f^dislocation.r * exp(-(dislocation.E + P * dislocation.V) / (dislocation.R * T))
    return Cdiff, Cdisl
end

function equal_mechanism_strain_rate(diffusion, dislocation, T, P, f, d)
    Cdiff, Cdisl = creep_prefactors(diffusion, dislocation, T, P, f, d)
    τ_eq = (Cdiff / Cdisl)^(1 / (dislocation.n - diffusion.n))
    ε̇_eq = 2 * Cdiff * τ_eq^diffusion.n
    return ε̇_eq, τ_eq
end

function equal_mechanism_surface(diffusion, dislocation;
        T_range = range(300.0, 2800.0, length = 121),
        P_range = range(0.1, 80, length = 121),
        f = 1.0,
        d = 1.0e-3)

    log10_ε̇II = zeros(length(T_range), length(P_range))
    τII = similar(log10_ε̇II)

    for (j, P_GPa) in pairs(P_range), (i, T_C) in pairs(T_range)
        T = T_C + 273.15
        P = P_GPa * 1.0e9
        ε̇_eq, τ_eq = equal_mechanism_strain_rate(diffusion, dislocation, T, P, f, d)

        log10_ε̇II[i, j] = log10(ε̇_eq)
        τII[i, j] = τ_eq
    end

    return (; T_range, P_range, log10_ε̇II, τII)
end

function plot_equal_mechanism_surface(map; d = 1.0e-3)
    fig = Figure(fontsize = 22, size = (1250, 850), backgroundcolor = :black)
    ax = Axis3(fig[1, 1],
        title = L"\mathrm{Diffusion-dislocation\ transition\ surface}",
        xlabel = L"T\ [^\circ\mathrm{C}]",
        ylabel = L"P\ [\mathrm{GPa}]",
        zlabel = L"\log_{10}\dot{\epsilon}_{II}\ [\mathrm{s}^{-1}]",
        backgroundcolor = :black,
        azimuth = 0.65π,
        elevation = 0.22π,
        protrusions = (70, 50, 80, 35),
        xlabeloffset = 55,
        ylabeloffset = 55,
        zlabeloffset = 70,
    )
    ax.titlecolor = :white
    ax.xlabelcolor = :white
    ax.ylabelcolor = :white
    ax.zlabelcolor = :white
    ax.xticklabelcolor = :white
    ax.yticklabelcolor = :white
    ax.zticklabelcolor = :white
    ax.xspinecolor_1 = :white
    ax.xspinecolor_2 = :white
    ax.xspinecolor_3 = :white
    ax.yspinecolor_1 = :white
    ax.yspinecolor_2 = :white
    ax.yspinecolor_3 = :white
    ax.zspinecolor_1 = :white
    ax.zspinecolor_2 = :white
    ax.zspinecolor_3 = :white
    ax.xgridcolor = (:white, 0.18)
    ax.ygridcolor = (:white, 0.18)
    ax.zgridcolor = (:white, 0.18)
    colsize!(fig.layout, 1, Relative(0.82))
    colgap!(fig.layout, 35)
    rowsize!(fig.layout, 1, Relative(0.88))

    sf = surface!(ax, map.T_range, map.P_range, map.log10_ε̇II;
        color = map.τII ./ 1.0e6,
        colormap = :lipari,
        shading = true,
        transparency = false,
    )

    cb = Colorbar(fig[1, 2], sf, label = L"\tau_{II}\ [\mathrm{MPa}]", width = 22, labelpadding = 16, ticklabelpad = 8)
    cb.labelcolor = :white
    cb.ticklabelcolor = :white
    cb.tickcolor = :white
    cb.spinewidth = 1
    cb.topspinecolor = :white
    cb.rightspinecolor = :white
    cb.leftspinecolor = :white
    cb.bottomspinecolor = :white
    Label(fig[2, 1:2],
        latexstring("\\mathrm{Surface\\ marks\\ dislocation\\ creep\\ fraction}=0.5;\\ \\mathrm{grain\\ size}=$(d * 1.0e3)\\ \\mathrm{mm}"),
        fontsize = 18,
        color = :white,
    )

    display(fig)
    return fig
end

diffusion, dislocation = let
    R = 8.314462618

    diffusion = DiffusionCreep(
        1,       # n, diffusion creep
        0.0,     # water fugacity exponent
        3.0,     # grain-size exponent
        3e-21,   # material parameter, SI-style example value
        160e3,   # activation energy [J mol^-1]
        8e-6,     # activation volume [m^3 mol^-1]
        R,
    )

    dislocation = DislocationCreep(
        3.0,     # n, dislocation creep
        0.0,     # water fugacity exponent
        7e-26,   # material parameter, SI-style example value
        190e3,   # activation energy [J mol^-1]
        8e-6,     # activation volume [m^3 mol^-1]
        R,
    )

    diffusion, dislocation
end

d = 1.0e-3
map = equal_mechanism_surface(diffusion, dislocation; d = d)
fig = plot_equal_mechanism_surface(map; d = d)
save("def_mechanism.png", fig)
