using Plots, Printf, LinearAlgebra, Statistics

function BraunBercoRic()
    # Main 
    H     = 1e3      # size
    Tbg   = 650+273 # temperature
    P     = 5e9
    Ebg   = 1e-15   # strain rate
    Vbg   = 2*H*Ebg # velocity
    dbg   = 1e-3    # grain size 
    G     = 6e10    # Shear modulus
    k     = 3
    rho   = 3300
    c     = 1000
    R     = 8.314
    # 0.1 Set material properties - same as in Bercovici et al 2018
    #----------------------------------------------
    # dislocation creep rheology
    n     = 3.5
    Adis  = 1.1*10^5*10^(-6*n) # prefactor in MPa and mum, Hirth & Kohlstedt 2003 --> converted to SI
    Qdis  = 530e3 # activation energy (gas constant included)
    Vadis = 14e-6

    # diffusion creep rheology
    m     = 3.0
    Adif  = 1.5e9.*10^(-6*m).*10^(-6) # Kohlstedt: 1.5e9, Hansen: 10^7.6, Karato & Wu: 8.7*10^(-15) 1/s
    Qdif  = 300e3
    Vadif = 6e-6

    # grain size evolution
    phi1    = 0.4 # 40 # pyroxene
    phi2    = 1-phi1
    eta_gs  = 3*phi1*phi2
    sigma   = 0.8 # half width of lognormal grain size distributions
    gamma   = 1
    c_gs    = 3/160 * exp(6*sigma^2)
    tmp     = phi1/sqrt(c_gs*(1-phi1)) +  phi2/sqrt(c_gs*(1-phi2)) # prefactor for conversion from grain curvature to grain size, should be approx pi/2 for phi1 = 0.4

    Ag      = 2*10^4 # Karato growth parameters --> converted to SI ?
    Qg      = 300e3
    Vag     = 6e-64
    p       = 2
    q       = 4 # this is from fitting Hiragas data with the 2-phase growth formulations (Bercovici&Ricard 2012)

    f_I         = 1e-4 # fraction of deformational work going into grain size reduction
    D_gse       = f_I/(gamma*eta_gs) # prefactor for grain size reduction as a combination of all other constant prefactors

    # temperature-dependent prefactors 
    # Bercovici et al.(2018) computed those prefactors taking into account a
    # pressure of 5 GPa and activation energies for both diffusion and
    # dislocation creep. Backwards engineering indicates that the following
    # activation energies were used: --> in the end, this equals a variation in
    # activation energies
    Hdis_eff = (Qdis+Vadis*P)
    Hdif_eff = (Qdif+Vadif*P)
    Hg_eff   = Hdif_eff

    # Space discretisation
    ncy    = 100
    dy     = H/ncy
    yv     = LinRange(0, H, ncy+1)
    yce    = LinRange(0-dy/2, H+dy/2, ncy+2)
    yc     = yce[2:end-1]
    Vxe    = collect(2*yce.*Ebg)
    Vx     = Vxe[2:end-1] 

    # Time discretisation
    nt     = 1000
    dtmax  = 1e11
    time   = 0

    # Initial arrays
    eta     = 1e22*ones(ncy+1)
    Q_dl    = zeros(ncy+1)
    Q_gss   = zeros(ncy+1)
    Q_dl_c  = zeros(ncy)
    Q_gss_c = zeros(ncy)
    eta_dl  = copy(eta)
    eta_gss = copy(eta)
    dred   = zeros(ncy+1)
    dgro   = zeros(ncy+1)
    A      = zeros(ncy+1)
    B      = zeros(ncy+1)
    GI     = zeros(ncy+1)
    G_gse  = zeros(ncy+1)
    # d      = dbg .* ( 1 .+ 0.01* cos.(2*π.*yv./H*8) .* exp.(-(yv.-H/2).^2/2/(H/2)^2) .* rand(ncy+1) )
    d      = dbg .* ( 1 .+ 0.07* cos.(2*π.*yv./H*8) .* exp.(-(yv.-H/2).^2/2/(H/2)^2) .* rand(ncy+1) )
    d0     = copy(d)
    d1     = copy(d)
    T      = Tbg*ones(ncy)
    Tv     = zeros(ncy+1)
    T0     = copy(T)
    Te     = Tbg*ones(ncy+2)
    qT     = zeros(ncy+1)
    Exy    = 1/2*diff(Vxe)/dy
    Eii    = copy(Exy)
    Txy    = 0*2*eta.*Exy
    Tii    = copy(Txy)
    Txy0   = copy(Txy)

    Fv     = zeros(ncy)
    Ft     = zeros(ncy)
    Fd     = zeros(ncy+1)


    dVdtau  = zeros(ncy)
    dTdtau  = zeros(ncy)
    dtauV   = zeros(ncy)
    etaVx   = zeros(ncy)
    etai    = copy(eta)
    eta1    = copy(eta)
    Q_dli   = copy(Q_dl)
    di      = copy(d)   

    # char stress
    tauc   = 1e8
    lc     = 1e-7
    tc     = 1/Ebg
    Ma     = 1e6*365*24*3600
    cma    = 365*24*3600*100

    # Pseudo-transient parameteres
    rel    = 0.1
    Reopt  = 0.5*pi
    cfl    = 6
    rhoi   = cfl*Reopt/ncy
    alpi   = H^2*c^2*rho^2/(4*pi^2*k) 
    Kappa  = k/rho/c
    eps    = 1e-9

    # Storage
    str         = zeros(nt)
    mean_visc   = zeros(nt)
    min_srate   = zeros(nt)
    mean_srate  = zeros(nt)
    max_srate   = zeros(nt)
    meca        = zeros(nt)
    mean_stress = zeros(nt)
    nFv0, nFd0, nFt0 = 1., 1., 1. 

    # Time loop
    anim = @animate for it=1:nt
        # Past
        d0   .= d
        Txy0 .= Txy
        T0   .= T
        
        @. Tv    = 0.5*(Te[1:end-1]  + Te[2:end])
        @. A     = Adis*exp(-Hdis_eff./R./Tv) # dislocation creep prefactor including temperature effects
        @. B     = (tmp)^(-m) * Adif*exp(-Hdif_eff./R./Tv) # diffusion creep prefactor
        @. GI    = (q/p *exp(-Hg_eff./R./Tv)* (Ag/250)) *(1e-6)^q # constant 
        @. G_gse = eta_gs.*GI./q # prefactor for grain growth as a combination of all other constant prefactors

        # Time step
        # maximum available deformational work to set the time step
    #     # Q_dl         = 2*1/2*abs(Exy).^(1/n-1).*A1.^(-1/n).*Exy.^2 # maximum total available deformational work
        @. Q_dl            = abs(Exy./A).^(1/n) .* Exy # maximum total available deformational work in dis creep
        @. Q_gss           = abs(Exy./B).*d.^m .* Exy
        @. dred            = -D_gse.*d.^2 .* (Q_dl + Q_gss)    
        dt           = 25*minimum(d)/maximum(abs.(dred))
        if dt>dtmax dt = dtmax end
        time         = time + dt
        @printf("Step %03d, dt = %02.2e \n", it, dt)
        # Iteration loop
        for iter=1:500000
            
            @. etai         = eta 
            @. Q_dli        = Q_dl
            @. di           = d   
            Vxe[2:end-1] = Vx                          # velocity + BC
            Vxe[1]       = 0 - Vx[1]                   # velocity + BC
            Vxe[end]     = 2*Vbg - Vx[end]             # velocity + BC
            Te[2:end-1]  = T                           # temperature + BC
            Te[1]        = T[1]                        # temperature + BC
            Te[end]      = T[end]                      # temperature + BC

            Exy          .= (1/2)*diff(Vxe)/dy            # strain rate
            @. Eii          = abs(Exy)                    # strain rate invariant
            @. Txy          = 2*eta.*(Exy + Txy0./(2*G*dt)) 
            @. Tii          = abs(Txy)                    # stress invariant using previous iteration guess
            @. eta_dl       = 1/2*Tii.^(1-n)/A            # from eq 2
            @. eta_gss      = 1/2*d.^m./B                 # from eq 3
            @. eta1         = 1 ./ (1 ./ eta_gss + 1 ./ eta_dl + 1 ./ (G*dt)) # eq 1          without factor
            @. eta          = rel.*eta1 + (1-rel)*etai    # relax non-linearity
            @. Q_dl         = Tii.^2 ./ 2 ./ eta_dl
            @. Q_gss        = Tii.^2 ./ 2 ./ eta_gss
            # Mechanical
            @. Txy          = 2*eta.*(Exy + Txy0./(2*G*dt))  # stress
            Fv             .= diff(Txy)/dy                # force balance
            # Thermal
            qT             .= -k*diff(Te)/dy
            @. Q_dl_c       = 0.5*( Q_dl[1:end-1]  + Q_dl[2:end])
            @. Q_gss_c      = 0.5*(Q_gss[1:end-1] + Q_gss[2:end])
            Ft          .= (1-f_I)*(Q_dl_c + Q_gss_c) - diff(qT)/dy - rho*c*(T-T0)/dt 
            # GSE
            @. dred         = -D_gse.*d.^2 * (Q_dl + Q_gss)      # gs reduction
            @. dgro         = G_gse./d.^(q-1)                    # gs growth
            
            @. d1           = dt*(dred + dgro) + d0
            @. d            = rel*d1 + (1-rel)*di
            @. Fd           = dt*(dred + dgro) + d0 - d 
            
            # Check residuals
            nFv          = norm(Fv)/length(Fv) / (tauc/lc)
            nFd          = norm(Fd)/length(Fd) / (lc)
            nFt          = norm(Ft)/length(Ft) / (tauc/c)
            if iter==1 nFv0 = nFv; nFd0 = nFd; nFt0 = nFt end
            if mod(iter,2000)==0 || iter==1
                @printf("It. %05d, F = %2.2e %2.2e - Fd = %2.2e %2.2e - Ft = %2.2e %2.2e\n", iter, nFv, nFv/nFv0, nFd, nFd/nFd0, nFt, nFt/nFt0)
                if  (nFv<eps || nFv/nFv0<eps)  && (nFd<eps || (nFd/nFd0)<eps) && (nFt<eps || (nFt/nFt0)<eps)             
                    break
                end
            end
            # Pseudo-transient solver
            @. etaVx  = 0.5*(eta[2:end] + eta[1:end-1])
            @. dtauV  = rhoi.*dy^2 ./etaVx ./2.1 .* cfl
            dtaut  = rhoi*dy^2/Kappa./2.1 .* cfl / 70e6
    # #         L=H
    # #         kT=k
    # #         
    # #         rhoc_i = 1
    # #         for i = 1:10
    # #                 dtaut     = rhoc_i*dy^2/Kappa./2.1 .* cfl / 1e35
    # #     
    # #         dtau = dtaut
    # #         rhoc_i = (L.^2.*c.*rho + 2*pi^2*dtau.*kT + 2*pi*sqrt(dtau.*kT.*(L.^2.*c.*rho - L.^2 + pi^2*dtau.*kT)))./L.^2
    # #         end

            dVdtau    = (1-rhoi) .* dVdtau + Fv
            dTdtau    = (1-rhoi) .* dTdtau + Ft
            Vx        = Vx + dtauV.*dVdtau
            T         = T  + dtaut.*dTdtau
        end

    # #     max(Q_dl)
    # #     fprintf( '#2.4e #2.4e\n\n', dred(1)*dt, (d(1)-d0(1)))
    # #     fprintf( '#2.4e #2.4e\n\n', d(1), d0(1) + dred(1)*dt)

        # Store
        str[it]         = time*Ebg/2
        mean_visc[it]   = mean(eta)
        mean_stress[it] = mean(Tii)
        min_srate[it]   = minimum(log10.(Eii))
        mean_srate[it]  = mean(log10.(Eii))
        max_srate[it]   = maximum(log10.(Eii))
        E_dl            = A.*Tii.^n
        E_gss           = B.*Tii.*d.^m
        meca[it]        = mean(E_dl ./ Eii)

        if mod(it,10)==0 || it==1
            eta_diss = 1 ./ (1 ./ eta_gss + 1 ./ eta_dl)
            p1 = plot( title = "Vx @ $(@sprintf("%1.6f", time/Ma)) Ma", xlabel = "Vx [cm/y]", ylabel = "y [m]" )
            p1 = plot!(Vx*cma, yc, label="")
            p2 = plot( title = "Eii @ $(@sprintf("%1.6f", time/Ma)) Ma", xlabel = "log10 Exy [1/s]", ylabel = "y [m]" , xlims=(-20, -10) )
            p2 = plot!(log10.(Eii), yv, label="")
            p3 = plot( title = "eta @ $(@sprintf("%1.6f", time/Ma)) Ma", xlabel = "log10 eta [Pa.s]", ylabel = "y [m]", xlims=(15, 25) )
            p3 = plot!(log10.(eta_diss), yv, label="" )
            p4 = plot( title = "d @ $(@sprintf("%1.6f", time/Ma)) Ma", xlabel = "log10 d [m]", ylabel = "y [m]" , xlims=(-6, -2) )
            p4 = plot!(log10.(d), yv, label="")
            display(plot(p1, p2, p3, p4))
        end
    end
    gif(anim, "GSE_Localization1D.gif", fps = 40)
end

BraunBercoRic()