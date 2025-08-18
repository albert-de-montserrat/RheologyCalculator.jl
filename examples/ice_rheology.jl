using Plots, Printf
include("RheologyDefinitions.jl")

function ice_rheology()
    verbose = true
    # SWITCHES
    series = 0 # see LocalIter_v3.ipynb for derivation of having eps_basal and eps_gbs in parallel
    elastic = 1
    gse = 1
    growth = 1
    reduct = 1
    diff = 1
    disl = 1
    bas_gbs = 1
    von_mises = 1
    # physics
    R = 8.314
    P = 1.0e7       # [Pa]             <--- should be the input from the code
    d_o = 1.0e-3      # d: grain size[m] <--- should be the input from the code
    G = 3.5e9     # Elastic modulus (E/1+nu/2 = 1e9/(1+0.3)/2)
    Coh = 1.0e7       # Von Mises yield
    fric = 0 * π / 180   # friction angel
    p0 = 7.0e-8      # p0: pressure heating coef[K/Pa]
    # GSE - From Behn et al. 2020
    lam_disl = 0.005
    lam_gbs = 0.005
    Qgg = 42.0e3
    Kgg = 1.02e-8
    gam = 0.065
    p = 6.03
    c = 3.0
    # numerics
    nt = 1         # number of time steps
    dt = 1.0e7       # physical time step value
    litmax = 100
    tol = 1.0e-15
    nT = 100
    nE = 101
    # init
    eta_el = G * dt
    d = d_o
    # 2D space
    Tvec = LinRange(220, 273, nT)
    Eiivec = 10.0 .^ LinRange(-15, -3, nE)
    Eii = 0 * Tvec .+ Eiivec'
    T = Tvec .+ 0 * Eiivec'

    # Initial dev. stress
    Txx = zeros(size(Eii))
    Tyy = zeros(size(Eii))
    Tzz = zeros(size(Eii))
    Txy = zeros(size(Eii))
    Tii = zeros(size(Eii))
    Eii_diff = zeros(size(Eii))
    Eii_basal = zeros(size(Eii))
    Eii_gbs = zeros(size(Eii))
    Eii_disl = zeros(size(Eii))
    Eii_bas_gbs = zeros(size(Eii))
    Eii_vm = zeros(size(Eii))
    mech = zeros(size(T))
    eta_ve = zeros(size(T))

    # Strain rate components and invariant
    exx = Eii ./ 2 .* ones(size(T))
    eyy = -Eii ./ 2 .* ones(size(T))
    ezz = 0 * ones(size(T))
    exy = 0 * ones(size(T))
    eii = sqrt.(1 / 2 * (exx .^ 2 + eyy .^ 2 + ezz .^ 2) + exy .^ 2)

    ## --- Goldsby
    T_h = T .+ p0 .* P  # T_h: melting point adjusted temperature[K]
    # Diffusional Rheology (Goldsby and Kohlstedt 2001, eq.4 + Table 6)
    m_diff = 2
    n_diff = 1
    F_diff = (2^((2 * n_diff - 1) / n_diff))^-1 # simple shear epxeriments (not sure 100#)
    Dv = 9.1e-4    # Dv: exponential prefactor[m^2/s]
    Vm = 1.97e-5   # Vm: molar volume
    A_diff = 42.0 .* Vm ./ (R .* T_h) .* Dv
    Q_diff = 59.4e3    # Q1: diffusional activation energy[J/mol]
    # Basal Rheology (Goldsby and Kohlstedt 2001, eq.3 + Table 5)
    n_basal = 2.4
    F_basal = (2^((n_basal - 1) / n_basal) * 3^((n_basal + 1) / (2 * n_basal)))^-1 # pure shear epxeriments
    A_basal = 5.5e7 * 10^(-6 * n_basal) # MPa to Pa conversion
    Q_basal = 60.0e3
    # GBS Rheology (Goldsby and Kohlstedt 2001, eq.3 + Table 5)
    m_gbs = 1.4
    n_gbs = 1.7
    F_gbs = (2^((n_gbs - 1) / n_gbs) * 3^((n_gbs + 1) / (2 * n_gbs)))^-1
    A1_gbs = 3.9e-3 * 10^(-6 * n_gbs) # A_gbs: preexponential constant[m^1.4/s*Pa^1.8] ######### ici l'unite c'est Pa^(-n_gbs).s^(-1).m^(m_gbs)
    Q1_gbs = 49.0e3                 # Q1_gbs: lower activation energy[J/mol]
    Q2_gbs = 197.0e3                # Q2_gbs: higher activation energy[J/mol]
    Tstar_gbs = 257                  # Tstar_gbs: activation threshold[K] (-18C)
    # Dislocation Rheology (Goldsby and Kohlstedt 2001, eq.3 + Table 5)
    n_disl = 4
    F_disl = (2^((n_disl - 1) / n_disl) * 3^((n_disl + 1) / (2 * n_disl)))^-1 # pure shear epxeriments
    A1_disl = 4.0e4 * 10.0^(-6 * n_disl) # A_disl: preexponential constant[1/s*Pa^4]  ######### ici l'unite c'est Pa^(-n_disl).s^(-1)
    Q1_disl = 64.0e3               # Q1_disl: lower activation energy[J/mol]
    Q2_disl = 220.0e3              # Q2_disl: higher activation energy[J/mol]
    Tstar_disl = 255                # Tstar_gbs: activation threshold[K] (-15C)
    # Temperature dependence of activation energy ONLY for GBS and DISL (melting point)
    Tstar_gbs_h = Tstar_gbs + p0 * P # Tstar_h:  adjusted activation threshold[K]
    Tstar_disl_h = Tstar_disl + p0 * P # Tstar_h:  adjusted activation threshold[K]
    Q_gbs_h = zeros(size(T))    # Q_gbs_h:  melting point adjusted activation energy[J/mol]
    Q_disl_h = zeros(size(T))    # Q_disl_h: melting point adjusted activation energy[J/mol]
    @. Q_gbs_h[(T_h - Tstar_gbs_h) < 0] = Q1_gbs
    @. Q_gbs_h[(T_h - Tstar_gbs_h) >= 0] = Q2_gbs
    @. Q_disl_h[(T_h - Tstar_disl_h) < 0] = Q1_disl
    @. Q_disl_h[(T_h - Tstar_disl_h) >= 0] = Q2_disl
    As_gbs = A1_gbs .* exp(-Q1_gbs ./ (R * Tstar_gbs_h))
    As_disl = A1_disl .* exp(-Q1_disl ./ (R * Tstar_disl_h))

    # Pre-computations
    B_diff = @. F_diff .* A_diff .^ (-1 / n_diff) .* exp(Q_diff ./ (n_diff .* R .* T_h))                                  # ---> CHECK Pa.s^(1/n)
    B_basal = @. F_basal .* A_basal .^ (-1 / n_basal) .* exp(Q_basal ./ (n_basal .* R .* T_h))                                  # ---> CHECK Pa.s^(1/n)
    B_gbs = @. F_gbs .* As_gbs .^ (-1 / n_gbs) .* exp((Q_gbs_h ./ (n_gbs .* R)) .* ((1.0 ./ T_h) - (1.0 ./ Tstar_gbs_h)))   # ---> CHECK Pa.s^(1/n)
    B_disl = @. F_disl .* As_disl .^ (-1 / n_disl) .* exp((Q_disl_h ./ (n_disl .* R)) .* ((1.0 ./ T_h) - (1.0 ./ Tstar_disl_h)))   # ---> CHECK Pa.s^(1/n)
    C_diff = @. (2.0 .* B_diff) .^ (-n_diff)    # ---> CHECK Pa^(-n).s^(-1)
    C_basal = @. (2.0 .* B_basal) .^ (-n_basal)
    C_gbs = @. (2.0 .* B_gbs) .^ (-n_gbs)
    C_disl = @. (2.0 .* B_disl) .^ (-n_disl)

    # Plasticity
    syield = Coh * cos(fric) + P * sin(fric)

    # V-E-P sans yield - pas de if
    # E = Ev + Ee + Ep
    # E = Ev + Ee + <Ep> (Ep si F>0) - discontinu

    for it in 1:1  #nt
        # Old solutions
        Txx_o = Txx
        Tyy_o = Tyy
        Tzz_o = Tzz
        Txy_o = Txy
        Tii_o = @. sqrt(1 / 2 * (Txx .^ 2 + Tyy .^ 2 + Tzz .^ 2) + Txy .^ 2)
        d_o = d

        # Viscoelastic strain rate components and invariant
        Exx = exx + elastic .* Txx_o ./ 2.0 ./ eta_el
        Eyy = eyy + elastic .* Tyy_o ./ 2.0 ./ eta_el
        Ezz = ezz + elastic .* Tzz_o ./ 2.0 ./ eta_el
        Exy = exy + elastic .* Txy_o ./ 2.0 ./ eta_el
        Eii = @. sqrt(1 / 2 * (Exx .^ 2 + Eyy .^ 2 + Ezz .^ 2) + Exy .^ 2)

        # Isolated viscosities
        eta_diff = B_diff .* Eii .^ (1.0 ./ n_diff - 1) .* d .^ (m_diff ./ n_diff)  # ---> CHECK Pa.s
        eta_basal = B_basal .* Eii .^ (1.0 ./ n_basal - 1)                       # ---> CHECK Pa.s
        eta_gbs = B_gbs .* Eii .^ (1.0 ./ n_gbs - 1) .* d .^ (m_gbs ./ n_gbs)    # ---> CHECK Pa.s
        eta_disl = B_disl .* Eii .^ (1.0 ./ n_disl - 1)                       # ---> CHECK Pa.s
        eta_vm = @. syield ./ 2.0 / eii                                   # ---> CHECK Pa.s

        # Define initial guesses
        if series == 0
            eta_ve = @. 1.0 / (diff ./ eta_diff + bas_gbs .* ((eta_basal .* eta_gbs) .^ (-1) ./ (eta_basal + eta_gbs)) + disl ./ eta_disl + elastic ./ eta_el + von_mises ./ eta_vm)
        elseif series == 1
            eta_ve = @. 1.0 / (diff ./ eta_diff + bas_gbs ./ eta_gbs + bas_gbs ./ eta_basal + disl ./ eta_disl + elastic ./ eta_el + von_mises ./ eta_vm)
        end
        d = d_o

        verbose && @printf("Starting local iterations:\n")
        re_vec = zeros(litmax)
        rd_vec = zeros(litmax)
        rese0, resd0 = 0.0, 0.0

        # Local iterations
        lit = 0
        while lit < litmax

            lit += 1

            # Function evaluation at current effective viscosity
            Tii .= 2 * eta_ve .* Eii

            # Strain rate invariants
            Eii_diff = @.     diff * C_diff .* Tii .^ n_diff .* d .^ (-m_diff)    # ---> CHECK s^(-1)
            Eii_basal = @.  bas_gbs * C_basal .* Tii .^ n_basal                   # ---> CHECK s^(-1)
            Eii_gbs = @.  bas_gbs * C_gbs .* Tii .^ n_gbs .* d .^ (-m_gbs)     # ---> CHECK s^(-1)
            Eii_disl = @.     disl * C_disl .* Tii .^ n_disl
            Eii_bas_gbs = @. bas_gbs * (1.0 / Eii_basal + 1.0 / Eii_gbs) .^ (-1)
            Eii_vm = @. von_mises * Tii ./ 2.0 / eta_vm

            if series == 0
                Eii_vis_0 = Eii_diff + Eii_bas_gbs + Eii_disl + Eii_vm
            elseif series == 1
                Eii_vis_0 = Eii_diff + Eii_basal + Eii_gbs + Eii_disl + Eii_vm
            end
            # Rheological residuals
            r_eta = Eii - Eii_vis_0 - elastic * Tii / 2 / eta_el

            # Update individual deviatoric stress components
            Txx = 2 * eta_ve .* Exx
            Tyy = 2 * eta_ve .* Eyy
            Tzz = 2 * eta_ve .* Ezz
            Txy = 2 * eta_ve .* Exy

            # Work that lead to grain size reduction
            Wxx_disl = Txx .* Eii_disl .* Txx ./ Tii
            Wyy_disl = Tyy .* Eii_disl .* Tyy ./ Tii
            Wzz_disl = Tzz .* Eii_disl .* Tzz ./ Tii
            Wxy_disl = Txy .* Eii_disl .* Txy ./ Tii
            Wxx_gbs = Txx .* Eii_gbs .* Txx ./ Tii
            Wyy_gbs = Tyy .* Eii_gbs .* Tyy ./ Tii
            Wzz_gbs = Tzz .* Eii_gbs .* Tzz ./ Tii
            Wxy_gbs = Txy .* Eii_gbs .* Txy ./ Tii
            Q_disl = Wxx_disl + Wyy_disl + Wzz_disl + 2 * Wxy_disl          # travail de deformation qui reduit la taille de grain
            Q_gbs = Wxx_gbs + Wyy_gbs + Wzz_gbs + 2 * Wxy_gbs           # travail de deformation qui reduit la taille de grain

            # Reduction rate
            dred = @. - d .^ 2 ./ c ./ gam .* (lam_disl .* Q_disl + lam_gbs .* Q_gbs) # d grain reduction eqn. (12)
            # Growth rate
            dgg = @.  p .^ (-1) .* d .^ (1 - p) .* Kgg .* exp(-Qgg ./ R ./ T)                # d grain growth eq. (4)
            # Grain size evolution residual
            r_d = @. (d - d_o) ./ dt - gse .* (growth .* dgg + reduct .* dred)         # r_d   (residual of grain size equation)

            # Residual check
            rese = maximum(abs.(r_eta[:]))
            resd = maximum(abs.(r_d[:]))
            if (lit == 1)
                rese0 = rese; resd0 = resd
            end
            verbose && @printf("It. %02d, re_rel = %2.2e re_abs = %2.2e, rd_rel = %2.2e rd_abs = %2.2e\n", lit, rese / rese0, rese, resd / resd0, resd)
            re_vec[lit] = rese / rese0
            rd_vec[lit] = resd / resd0
            if (rese / rese0 < tol  || rese < tol) && (resd / resd0 < tol  || resd < tol)
                break
            end

            # Exact derivative
            if series == 0
                dr_eta_deta = @. -Eii .* elastic ./ eta_el - Eii .* von_mises ./ eta_vm - Eii_bas_gbs .* (Tii .^ (-n_gbs) .* d .^ m_gbs .* n_gbs ./ (C_gbs .* eta_ve) + Tii .^ (-n_basal) .* n_basal ./ (C_basal .* eta_ve)) ./ (Tii .^ (-n_gbs) .* d .^ m_gbs ./ C_gbs + Tii .^ (-n_basal) ./ C_basal) - Eii_diff .* n_diff ./ eta_ve - Eii_disl .* n_disl ./ eta_ve
                dr_eta_dd = @. Eii_diff .* m_diff ./ d + Eii_bas_gbs .* Tii .^ (-n_gbs) .* d .^ m_gbs .* m_gbs ./ (C_gbs .* d .* (Tii .^ (-n_gbs) .* d .^ m_gbs ./ C_gbs + Tii .^ (-n_basal) ./ C_basal))
                dr_d_deta = @. 2 * d .^ 2.0 * gse .* reduct .* (Eii_disl .* lam_disl .* (Exx .^ 2.0 * n_disl + Exx .^ 2 + 2 * Exy .^ 2.0 * n_disl + 2 * Exy .^ 2 + Eyy .^ 2.0 * n_disl + Eyy .^ 2 + Ezz .^ 2.0 * n_disl + Ezz .^ 2) + Eii_gbs .* lam_gbs .* (Exx .^ 2.0 * n_gbs + Exx .^ 2 + 2 * Exy .^ 2.0 * n_gbs + 2 * Exy .^ 2 + Eyy .^ 2.0 * n_gbs + Eyy .^ 2 + Ezz .^ 2.0 * n_gbs + Ezz .^ 2)) ./ (Eii .* c .* gam)
                dr_d_dd = @. (Eii .* c .* d .* gam + dt .* gse .* (Eii .* c .* dgg .* gam .* growth .* (p - 1) + 2 * d .^ 2.0 * eta_ve .* reduct .* (2 * Eii_disl .* lam_disl - Eii_gbs .* lam_gbs .* m_gbs + 2 * Eii_gbs .* lam_gbs) .* (Exx .^ 2 + 2 * Exy .^ 2 + Eyy .^ 2 + Ezz .^ 2))) ./ (Eii .* c .* d .* dt .* gam)
            elseif series == 1
                dr_eta_deta = @. -C_basal .* Tii .^ n_basal .* n_basal ./ eta_ve - C_gbs .* Tii .^ n_gbs .* d .^ (-m_gbs) .* n_gbs ./ eta_ve - Eii .* elastic ./ eta_el - Eii .* von_mises ./ eta_vm - Eii_diff .* n_diff ./ eta_ve - Eii_disl .* n_disl ./ eta_ve
                dr_eta_dd = @. C_gbs .* Tii .^ n_gbs .* d .^ (-m_gbs) .* m_gbs ./ d + Eii_diff .* m_diff ./ d
                dr_d_deta = @. 2 * d .^ 2.0 * gse .* reduct .* (Eii_disl .* lam_disl .* (Exx .^ 2.0 * n_disl + Exx .^ 2 + 2 * Exy .^ 2.0 * n_disl + 2 * Exy .^ 2 + Eyy .^ 2.0 * n_disl + Eyy .^ 2 + Ezz .^ 2.0 * n_disl + Ezz .^ 2) + Eii_gbs .* lam_gbs .* (Exx .^ 2.0 * n_gbs + Exx .^ 2 + 2 * Exy .^ 2.0 * n_gbs + 2 * Exy .^ 2 + Eyy .^ 2.0 * n_gbs + Eyy .^ 2 + Ezz .^ 2.0 * n_gbs + Ezz .^ 2)) ./ (Eii .* c .* gam)
                dr_d_dd = @. (Eii .* c .* d .* gam + dt .* gse .* (Eii .* c .* dgg .* gam .* growth .* (p - 1) + 2 * d .^ 2.0 * eta_ve .* reduct .* (2 * Eii_disl .* lam_disl - Eii_gbs .* lam_gbs .* m_gbs + 2 * Eii_gbs .* lam_gbs) .* (Exx .^ 2 + 2 * Exy .^ 2 + Eyy .^ 2 + Ezz .^ 2))) ./ (Eii .* c .* d .* dt .* gam)
            end
            # invert 2x2 Jacobian matrix
            det = @.  dr_eta_deta .* dr_d_dd - dr_eta_dd .* dr_d_deta
            D11 = @.  1.0 / det .* dr_d_dd
            D12 = @. -1.0 / det .* dr_eta_dd
            D21 = @. -1.0 / det .* dr_d_deta
            D22 = @.  1.0 / det .* dr_eta_deta
            deta = -(D11 .* r_eta + D12 .* r_d)                 # inverse times rhs r_eta is deta and r_d is fd
            dd = -(D21 .* r_eta + D22 .* r_d)

            # Newton update
            eta_ve = @. eta_ve + deta
            d = @. d + dd
        end #lit

        p1 = Plots.plot()
        p1 = Plots.plot!(1:lit, log10.(re_vec[1:lit]))
        p1 = Plots.plot!(1:lit, log10.(rd_vec[1:lit]))
        Plots.display(p1)

        # Dominant mechanisms
        Eii_el = elastic * (Tii - Tii_o) / 2 / eta_el
        if series == 1
            @. mech[Eii_diff > Eii_basal && Eii_diff > Eii_gbs   && Eii_diff > Eii_disl && Eii_diff > Eii_vm] = 1
            @. mech[Eii_basal > Eii_diff  && Eii_basal > Eii_gbs   && Eii_basal > Eii_disl && Eii_basal > Eii_vm] = 2
            @. mech[Eii_gbs > Eii_diff  && Eii_gbs > Eii_basal && Eii_gbs > Eii_disl && Eii_gbs > Eii_vm] = 3
            @. mech[Eii_disl > Eii_diff  && Eii_disl > Eii_basal && Eii_disl > Eii_gbs  && Eii_disl > Eii_vm] = 4
            @. mech[Eii_vm > Eii_diff  && Eii_vm > Eii_basal && Eii_vm > Eii_gbs  && Eii_vm > Eii_disl] = 5
        else
            eps_tmp = (1.0 ./ Eii_basal .+ 1.0 ./ Eii_gbs) .^ (-1)
            if diff == 1
                @. mech[Eii_diff > eps_tmp   && Eii_diff > Eii_disl && Eii_diff > Eii_el  && Eii_diff > Eii_el] = 1
            end
            if disl == 1
                @. mech[Eii_disl > Eii_diff  && Eii_disl > eps_tmp  && Eii_disl > Eii_el  && Eii_disl > Eii_el] = 4
            end
            if bas_gbs == 1
                @. mech[eps_tmp > Eii_diff  && eps_tmp > Eii_disl && eps_tmp > Eii_el  && eps_tmp > Eii_el] = 2
            end
            if von_mises == 1
                @. mech[Eii_vm > Eii_diff  && Eii_vm > Eii_disl && Eii_vm > eps_tmp && Eii_vm > Eii_el] = 5
            end
            @. mech[mech == 2 && Eii_gbs < Eii_basal] = 3
        end

    end

    p2 = Plots.heatmap(Tvec, log10.(Eiivec), log10.(d)', title = "log10 d", xlabel = "T [K]", ylabel = "log10 ε̇ [1/s]")
    p3 = Plots.heatmap(Tvec, log10.(Eiivec), log10.(Tii)', title = "log10 τII", xlabel = "T [K]", ylabel = "log10 ε̇ [1/s]")
    p4 = Plots.heatmap(Tvec, log10.(Eiivec), mech', title = "mechanism", xlabel = "T [K]", ylabel = "log10 ε̇ [1/s]")
    p5 = Plots.heatmap(Tvec, log10.(Eiivec), log10.(eta_ve)', title = "log10 η", xlabel = "T [K]", ylabel = "log10 ε̇ [1/s]")

    return Plots.display(Plots.plot(p2, p3, p4, p5))

end

ice_rheology()
