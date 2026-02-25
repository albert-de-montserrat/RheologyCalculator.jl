using GLMakie
using MathTeXEngine
Makie.update_theme!( fonts = (regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))

function PlotYieldFunction()

    τ = LinRange(0, 2e8, 1000)
    P = LinRange(-1e8, 5e8, 1500)

    fig = Figure(fontsize = 30, size = (800, 600) .* 1)
    ax1 = Axis(fig[1, 1], xlabel = L"$P$ [MPa]", ylabel = L"$\tau$ [MPa]", title=L"$$Yield functions", aspect=DataAspect())

    # Drucker-Prager
    C = 5e7
    ϕ = 35.0
    F = @. τ' - P*sind(ϕ) - C*cosd(ϕ)
    contour!(ax1, P./1e6, τ./1e6, F, label="Drucker-Prager", levels=[0], color=:black)

    # Fairhurst
    Pc = 300e6+1.4e8
    Pt = 10e6+1.4e8
    Pm = (Pt + Pc)/2
    qm = 2e8
    F = @. (τ'./qm).^2 + ((P - Pm)/(Pc - Pt))^2 - 1.0
    contour!(ax1, P./1e6, τ./1e6, F, label="Fairhurst", levels=[0], color=:blue)

    # Hyperbolic
    σt = 10e6
    F  = @. sqrt((τ').^2 + (C * cosd(ϕ) - σt*sind(ϕ))^2) - (P * sind(ϕ) + C * cosd(ϕ))
    contour!(ax1, P./1e6, τ./1e6, F, label="Hyperbolic", levels=[0], color=:red)

    # MCC
    Pt = -1e7
    M  = 0.9
    r  = 2e8
    β  = 0.1
    b  = zeros(size(F))
    b[P .+ 0.0.*τ' .<  Pt + r ] .= 1
    b[P .+ 0.0.*τ' .>= Pt + r ] .= β

    F  = @. 1/b *(P - Pt - r)^2  + (τ')^2 / M^2 - r^2 
    contour!(ax1, P./1e6, τ./1e6, F, label="Mod. Cam-Clay", levels=[0], color=:green)
    

    # M   = 0.8           # Critical state slope
    # p_c = 400e6#-10e6           # Preconsolidation pressure (Pa)
    # a   = 400e6           # Transition pressure (Pa)
    # b  = zeros(size(F))
    # β  = 1#0.1
    # b[P .+ 0.0.*τ' .<  a/2 ] .= 1
    # b[P .+ 0.0.*τ' .>=  a/2 ] .= β
    # F  = @. (τ' / (M))^2 + 1/b^2*(P - a)*(P - p_c)
    # contour!(ax1, P./1e6, τ./1e6, F, label="Mod. Cam-Clay", levels=[0], color=:yellow)


    # σt = 10e6
    # Pc = 4e8
    # M  =  0.7
    # F  = @. (τ').^2 / M^2 + P * (P - Pc)
    # contour!(ax1, P./1e6, τ./1e6, F, label="Mod. Cam-Clay", levels=[0], color=:green)
    

    # Hoek-Brown
    Pt = 4e8
    σt = 10e6
    M  = 0.9
    a  = 2e8
    b  = zeros(size(F))
   

    σCI = 1e8
    mi  = 10
    GSI = 50
    D   = 0.

    mb = mi*exp((GSI-100)/(28-14*D))
    s  = exp((GSI-100)/(9-3*D))
    a  = 1/2+1/6*(exp(-GSI/15) - exp(-20/3))

    F  = @. (τ') - σCI*abs(mb*(P - 1/3*(τ'))/σCI + s)^a
    contour!(ax1, P./1e6, τ./1e6, F, label="Hoek-Brown", levels=[0], color=:purple)
    

    axislegend(ax1, position = :lt, labelsize=12)
    display(fig)

end

PlotYieldFunction()