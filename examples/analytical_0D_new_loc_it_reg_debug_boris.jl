using Plots
using Roots # root-finding algorith
using ForwardDiff
using LinearAlgebra,  LaTeXStrings

struct MatProp  # will hold material properties
   #linear viscosity
   eta
   # non-linear viscosity
   Bn
   n
   # elasticity
   G
   Kb
   # plasticity
   ch
   fr
   dil
   Pt
   # regularization
   eta_min
   eta_vp
   # temporary
   ty  
end

function DevRes(Eps_II,Mat,dt,lambda_dot,dQdτ,Tau_II,TauII_0)

    if Mat.G == 0
        A_e    = 0.0;
        Eps_el = 0.0;
    else
        A_e     =   1.0/(2.0*Mat.G*dt);
        Eps_el  =  (Tau_II-TauII_0)*A_e;
    end
    if Mat.eta == 0
        A_l   = 0.0;
        Eps_l = 0.0;
    else
        A_l     =   1.0/(2.0*Mat.eta);
        Eps_l   =   Tau_II.*A_l;                    # Linear viscous strainrate  
    end
    if Mat.Bn == 0
        A_n   = 0.0;
        n     = Mat.n;
        Eps_n = 0.0;
    else
        A_n     =   Mat.Bn;
        n       =   Mat.n;
        Eps_n   =   A_n*Tau_II.^n;                  # Powerlaw viscosity strain rate
    end

    Eps_pl = lambda_dot*dQdτ

    #Eps_II  =   Eps_II+1e-9;

    r = Eps_II - Eps_el - Eps_l - Eps_n - Eps_pl

    return r, A_e, A_l, A_n

end

function BulkRes(Eps_vol,P_tr,P_o,Mat,dt,lambda_dot,dQdP)

    Eps_el = - (P_tr-P_o)/Mat.Kb/dt
    Eps_pl = lambda_dot*dQdP

    r = Eps_vol - Eps_el - Eps_pl

    #@show Eps_vol,Eps_el,Eps_pl,P_tr

    return r

end

function GetRes(xx,Eps_II,Eps_vol,P_o,Mat,dt,F,Tau_o_II,k, kf, c, a, b, pd, py, pf, R, flag, plast)

    Tau_II     = xx[1]
    P          = xx[2]
    lambda_dot = xx[3]

    F,dFdτ,dFdP,dQdτ,dQdP,dQdτdτ,dQdτdP,dQdPdτ,dQdPdP,st,flag = getPlastParam(k, kf, c, a, b, pd, py, pf, R, P, Tau_II,flag,plast)
    #r = zeros(3)

    r1, A_e, A_l, A_n = DevRes(Eps_II,Mat,dt,lambda_dot,dQdτ,Tau_II,Tau_o_II)
    r2 = BulkRes(Eps_vol,P,P_o,Mat,dt,lambda_dot,dQdP)
    #r3 = F - Mat.eta_vp * lambda_dot
    r3 = F - lambda_dot
    

    r = [r1,r2,r3]

    return r

end


# Residual vector as ordered in GP:
function GetRes1(xx,Eps_II,Eps_vol,P_o,Mat,dt,F,Tau_o_II,k, kf, c, a, b, pd, py, pf, R, flag, plast)

    Tau_II     = xx[1]
    lambda_dot = xx[2]
    P          = xx[3]
    
    F,dFdτ,dFdP,dQdτ,dQdP,dQdτdτ,dQdτdP,dQdPdτ,dQdPdP,st,flag = getPlastParam(k, kf, c, a, b, pd, py, pf, R, P, Tau_II,flag,plast)
    #r = zeros(3)

    r1, A_e, A_l, A_n = DevRes(Eps_II,Mat,dt,lambda_dot,dQdτ,Tau_II,Tau_o_II)
    r2 = F - Mat.eta_vp * lambda_dot
    #r2 = F - lambda_dot
    
    r3 = BulkRes(Eps_vol,P,P_o,Mat,dt,lambda_dot,dQdP)

   # @show F
    r = [r1,r2,r3]

    return r

end

function get_yield_param(Mat)

    # DP and flow potential parameters
    k   = sin(Mat.fr)
    kf  = sin(Mat.dil)
    c   = cos(Mat.fr)*Mat.ch
    
    # cap center, radius and scaling
    a    = sqrt(1 + k*k)
    cosa = 1/a
    sina = k/a
    py   = (Mat.Pt + c*cosa)/(1 - sina)
    R    = py - Mat.Pt

    # cap-DP delimiter
    pd = py - R*sina
    
    # flow potential center and scaling
    pf = pd + kf*(c + k*pd)
    b  = sqrt(1 + kf*kf)
    
    # return parameters
    return  k, kf, c, a, b, pd, py, pf, R

end

function getPlastParam(k, kf, c, a, b, pd, py, pf, R, P, Tau_II,flag,plast)

    #stabilization for local jacobian
    st = 0.

    # delimiter point
    τd = c + k*pd
    
    # yield function 
    if Tau_II > (py - P)*τd/(py - pd)
        
        F    =  Tau_II - k*P - c
        dFdτ =  1.
        dFdP = -k
        
    else
        
        Ry   = sqrt(Tau_II^2 + (P - py)^2)
        F    = a*(Ry - R)
        dFdτ = a/Ry*Tau_II
        dFdP = a/Ry*(P-py)

    end

    # flow potential 
    if Tau_II > (pf - P)*τd/(pf - pd)

        dQdτ   = 1. / 2.
        dQdP   = kf
        dQdτdτ = 0.
        dQdτdP = 0.
        dQdPdτ = 0.
        dQdPdP = 0.

    else
        
        Rf     =  sqrt(Tau_II^2 + (P - pf)^2)
        dQdτ   =  1. / 2. *b/Rf*Tau_II  
        dQdP   = -b/Rf*(P-pf)
        dQdτdτ =  1. / 2. *b*(Rf^2-Tau_II^2)/Rf^3
        dQdτdP = -b*Tau_II*(P - pf)/2. /Rf^3
        dQdPdτ =  b*Tau_II*(P - pf)/Rf^3
        dQdPdP = -b * (Rf^2 - (P - pf)^2)/Rf^3

    end

    if F < 0 && flag || plast == 0
        F      = 0.
        dFdτ   = 0.
        dFdP   = 0.
        dQdτ   = 0.
        dQdP   = 0.
        dQdτdτ = 0.
        dQdτdP = 0.
        dQdPdτ = 0.
        dQdPdP = 0.
        st     = 1.
    else
        flag = false
    end

    return F,dFdτ,dFdP,dQdτ,dQdP,dQdτdτ,dQdτdP,dQdPdτ,dQdPdP,st,flag

end

function GetJacLocalIt(xx,A_e, A_l, A_n,dFdτ,dFdP,dQdτ,dQdP,dQdτdτ,dQdτdP,dQdPdτ,dQdPdP,Mat,dt,st)
    # Note: this is reordered w.r.t. the paper to have [τ, λ, P]
    Tau_II     = xx[1]
    lambda_dot = xx[2]

    J = zeros(3,3)

    J[1, 1] = -lambda_dot*dQdτdτ -A_l -A_e  - Mat.n*A_n*Tau_II^(Mat.n-1.) 
    J[1, 2] = -dQdτ
    J[1, 3] = -lambda_dot*dQdτdP
    
    J[2, 1] =  dFdτ
    J[2, 2] =  -Mat.eta_vp
    J[2, 3] =  dFdP
    
    J[3, 1] = -lambda_dot*dQdPdτ
    J[3, 2] = -dQdP
    J[3, 3] = -lambda_dot*dQdPdP + 1. /Mat.Kb/dt
    
    
    return J

end

#--------------------------------------------------------------------------
function LocalIterations(Tau_o_II, Eps_II, Mat, dt, P, P_o, Eps_vol)
#
# Perform local iterations to determine Tau_II; this can be done with fzero 
# (root finding method) or with Newton-Rhapson iterations


    # initial guess
    Tau_II = Tau_o_II

    r          = zeros(3)
    xx         = zeros(3)
    rn         = zeros(3)
    rnorm      = 1
    A_e        = 0
    tol        = 1e-5
    maxit      = 1e3
    lambda_dot = 0
    flag       = true
    check      = 0
    plast      = 0
    k, kf, c, a, b, pd, py, pf, R = get_yield_param(Mat)

    while check == 0

        for i = 1:maxit

            F,dFdτ,dFdP,dQdτ,dQdP,dQdτdτ,dQdτdP,dQdPdτ,dQdPdP,st,flag = getPlastParam(k, kf, c, a, b, pd, py, pf, R, P, Tau_II,flag,plast)

            #create solution vector 

            xx = [Tau_II,lambda_dot, P]

            # get residual
            #@show F, flag, plast
            r = GetRes1(xx,Eps_II,Eps_vol,P_o,Mat,dt,F,Tau_o_II,k, kf, c, a, b, pd, py, pf, R, flag,plast)

            #@show r

            # get normalized residual
            rn[1] = abs(r[1]/Eps_II)
            rn[2] = abs(r[2]/c)
            rn[3] = abs(r[3]/Eps_vol)
            
            if Eps_vol == 0
                rn[3] = r[3]
            end

            rnorm = max(rn...)
           
            #@show r,rn,rnorm,Eps_II,Eps_vol,c

            if rnorm <= tol
                break
            end
#
            # Analytical jacobian:
            #_, A_e, A_l, A_n = DevRes(Eps_II,Mat,dt,lambda_dot,dQdτ,Tau_II,Tau_o_II)
            #J = GetJacLocalIt(xx,A_e, A_l, A_n,dFdτ,dFdP,dQdτ,dQdP,dQdτdτ,dQdτdP,dQdPdτ,dQdPdP,Mat,dt,st)

           # @show J, dFdτ,dFdP,dQdτ,dQdP,dQdτdτ,dQdτdP,dQdPdτ,dQdPdP

            #Check jacobian with AD
            f_res = (x -> GetRes1(x,Eps_II,Eps_vol,P_o,Mat,dt,F,Tau_o_II, k, kf, c, a, b, pd, py, pf, R, flag,plast))
            #@show f_res
            J   = ForwardDiff.jacobian(f_res,xx)
            
            #@show J
#            error("stop")
            #Check jacobian with Finite difference

            #jac = zeros(3,3)
            #eh  = zeros(3)
            #h   = 1e-8
#   
            #for j = 1:3
#   
            #    eh[j]    = h
            #    jac[:,j] = (GetRes(xx + eh,Eps_II,Eps_vol,P_o,Mat,dt,F,Tau_o_II,k, kf, c, a, b, pd, py, pf, R, flag,plast) - r)/h
            #    eh[j]    = 0.
            #end

            #if sum(abs.(jac[3,:])) < 1e-10
            #    jac[3,3] = 1.
            #end

            #@show jac

            # solve for update
            x  = J\-r

            #@show i, x
            # update variables
            Tau_II     += x[1]
            lambda_dot += x[2]
            P          += x[3]

        end

        F,dFdτ,dFdP,dQdτ,dQdP,dQdτdτ,dQdτdP,dQdPdτ,dQdPdP,st,flag = getPlastParam(k, kf, c, a, b, pd, py, pf, R, P, Tau_II,flag,1)

        if F < 1e-8 || plast == 1
            check = 1
            #@show plast, F
        else
            plast = 1
            #if its == 11
            #print("plasticity on")
            #end
        end

    end

    #@show flag

    if rnorm > tol
        div = 1
        @show div
    else
        div = 0
        #@show div
    end


    Eps_eff =   Eps_II + 1e-9 + Tau_o_II*A_e;

    Eta_eff  =   Tau_II/2/Eps_eff;

return Tau_II, Eta_eff, lambda_dot, P
end


function ComputeStressUpdate(Eps, P_tr, Tau_o, Mat, dt, P_o)
    # This computes the stress update, given old dev. stresses and new/old pressure
    
    # Perform local (scalar) iterations to obtain the converged stresses Tau_II and the effective viscosity:
    Iden        =   Matrix(LinearAlgebra.I, 4, 4)*1.0;
    m           =   [1;1;0;0;];
    Id          =   Iden-(1.0/2.0)*m*m';
    
    Tau_II_o    =   sqrt(sum(0.5*(Tau_o.*Tau_o)));  # old deviatoric stress invariant
    Eps_vol     =   sum(m.*Eps);                    # Volumetric strainrate (2D)
    d           =   Id*Eps;                         # Deviatoric strain rate
   
    d_II        =   sqrt(sum(0.5*(d  .*d  )));      # 2nd invariant dev. strainrate
    Eps_II      =   sqrt(sum(0.5*(Eps  .*Eps  )));  # 2nd invariant total strainrate 
    
    # Perform local iterations for nonlinear shear rheology
    Tau_II, Eta_eff, lambda_dot, P =   LocalIterations(Tau_II_o, d_II, Mat, dt, P_tr, P_o, Eps_vol);

    # Compute stress update  
    d_star      =   d + Tau_o/(2*Mat.G*dt);
    d_II_star   =   sqrt(sum(0.5*(d_star  .*d_star  )))
    Tau      =   2*Eta_eff .* d_star;            # new deviatoric (trial) stress

    return Tau, Eta_eff, lambda_dot, P

end

function AltComputeStressUpdate(d_II, Eps_vol, P_tr, Tau_II_o, Mat, dt, P_o)
    # This computes the stress update, given old dev. stresses and new/old pressure
    # Perform local iterations for nonlinear shear rheology
    Tau_II, Eta_eff, lambda_dot, P =   LocalIterations(Tau_II_o, d_II, Mat, dt, P_tr, P_o, Eps_vol);

    return Tau_II, Eta_eff, lambda_dot, P

end

function PlottingYield(Mat,P_Y)

    # Yield function

    # Cap function
    a       =   sqrt(1+sin(Mat.fr)*sin(Mat.fr))
    cosa    =   1/a
    sina    =   sin(Mat.fr)/a
    Py      =   (Mat.Pt+Mat.ch*cos(Mat.fr)*cosa)/(1-sina)
    R       =   Py-Mat.Pt
    Pj      =   Py - sin(atan(sin(Mat.fr)))*R

    P_cap   =   P_Y[P_Y .<= Pj]

    Tau_y_C2=   .-(P_cap .- Py).^2 .+ R^2
    Tau_y_C2[sign.(Tau_y_C2) .== -1] .= 0
    Tau_y_C =   sqrt.(Tau_y_C2)


    # Drucker-Prager
    P_DP    =   P_Y[P_Y .>= Pj]
    Tau_y_DP =   Mat.ch .* cos(Mat.fr) .+ P_DP .* sin(Mat.fr);

    Tauj    = Mat.ch .* cos(Mat.fr) .+ Pj .* sin(Mat.fr);

    Tau_y    =   vcat(Tau_y_C,Tau_y_DP);

    return Tau_y,Pj,Tauj,Py

end

function MainCode(d_II,Eps_vol,dt)
    # Material parameters:
    G             =   1e10
    Kb            =   2e11
    eta           =   1e20        # linear viscosity
    Bn            =   0;          # powerlaw prefactor
    n             =   2;          # powerlaw exponent
    fr            =   30/180*pi;  # friction angle
    dil           =   10/180*pi;   # dilation angle
    ch            =   1.0e6 #*1e10;    # cohesion
    Pt            =   -0.5e6         # tensile strength
    eta_min       =   1e-20
    eta_vp        =   1e-19
    ty            =   "DP"

    Mat     =   MatProp(eta,Bn,n,G,Kb,ch,fr,dil,Pt,eta_min,eta_vp,ty)

    # Applied background deviatoric strainrates

    P_o      =    0.0;
    Tau_II_o =    0.0

    # Numerics
    #dt          =   t_Max/30;
    nt          =   10;
    
    time         =   0;
    time_vec     =   zeros(nt+1);
    P_vec        =   zeros(nt+1); P_vec[1]=P_o;
    T2nd_vec     =   zeros(nt+1); T2nd_vec[1]=Tau_II_o;
    for itime=2:nt+1
    
        P_tr    =   P_o - dt*Mat.Kb*Eps_vol;     # trial P update (no plasticity)
        
        # Update deviatoric stresses, theta (and pressure in case of dilatant plasticity)
        Tau_II, Eta_eff, lambda_dot,  P = AltComputeStressUpdate(d_II, Eps_vol, P_tr, Tau_II_o, Mat, dt, P_o)
        #AltComputeStressUpdate(d_II, Eps_vol, P_tr, Tau_II_o, Mat, dt, P_o)
        #@show Eta_eff
       
        # Extract P and Tau (for visualization)
    
        time                =   time + dt
        time_vec[itime]     =   time
        T2nd_vec[itime]   =   Tau_II
        P_vec[itime]        =   P
        
        Tau_II_o   =   Tau_II
        P_o     =   P

        @show itime
        
    end


    return time_vec, T2nd_vec, P_vec, Mat, Pt
end

#CharDim  = GEO_units(length=10km, temperature=1000C, stress=10MPa, viscosity=1e20Pas)
SecYear  =    3600*24*365.25
d_II     =    0.0
Eps_vol  =    7e-15;
dt       =    2*SecYear #nondimensionalize(5e-5Myrs,CharDim);

time_vec1, T2nd_vec1, P_vec1, Mat, Pt = MainCode(d_II,Eps_vol,dt)

#=
d_II     =    nondimensionalize(7e-14s^-1,CharDim)
Eps_vol  =    0.0
dt       =   nondimensionalize(2.5e-6Myrs,CharDim)

time_vec2, T2nd_vec2, P_vec2 = MainCode(d_II,Eps_vol,dt,CharDim)

d_II     =    nondimensionalize(7e-14s^-1,CharDim)
Eps_vol  =    nondimensionalize(7e-15s^-1,CharDim)
dt       =   nondimensionalize(1.5e-6Myrs,CharDim)

time_vec3, T2nd_vec3, P_vec3 = MainCode(d_II,Eps_vol,dt,CharDim)  

P_Y = Pt:1e-3:0.2
Tau_y,Pj,Tauj,Py   = PlottingYield(Mat,P_Y)

time_vec1 = getfield.(dimensionalize(time_vec1,yr,CharDim),1)
T2nd_vec1 = getfield.(dimensionalize(T2nd_vec1,MPa,CharDim),1)
P_vec1    = getfield.(dimensionalize(P_vec1,MPa,CharDim),1)
time_vec2 = getfield.(dimensionalize(time_vec2,yr,CharDim),1)
T2nd_vec2 = getfield.(dimensionalize(T2nd_vec2,MPa,CharDim),1)
P_vec2    = getfield.(dimensionalize(P_vec2,MPa,CharDim),1)
time_vec3 = getfield.(dimensionalize(time_vec3,yr,CharDim),1)
T2nd_vec3 = getfield.(dimensionalize(T2nd_vec3,MPa,CharDim),1)
P_vec3    = getfield.(dimensionalize(P_vec3,MPa,CharDim),1)
P_Y       = getfield.(dimensionalize(P_Y,MPa,CharDim),1)
Tau_y     = getfield.(dimensionalize(Tau_y,MPa,CharDim),1)
Pj        = getfield.(dimensionalize(Pj,MPa,CharDim),1)
Py        = getfield.(dimensionalize(Py,MPa,CharDim),1)
Tauj      = getfield.(dimensionalize(Tauj,MPa,CharDim),1)


=#

# plot
pl = Plots.plot(time_vec1/SecYear, T2nd_vec1/1e6, label=L"τ_\mathrm{II}",lw = 1.5, xlabel="Time [yrs]", ylabel="Stress [MPa]", linecolor =:steelblue1, framestyle = :box);
pl = Plots.plot!(time_vec1/SecYear, P_vec1/1e6, label=L"P",lw = 2,legend=:bottomleft, linecolor =:indianred1, tickfontsize = 13, labelfontsize	= 15, legendfontsize = 12);
display(pl)
savefig("Time_evolution1.svg")


#=
pl1 = plot(P_vec1, T2nd_vec1, markershape=:circle, label="Pure extension",lw = 1.5, xlabel=L"P \; \mathrm{[MPa]}", ylabel=L"τ_\mathrm{II} \; \mathrm{[MPa]}", linecolor =:gold1, framestyle = :box, markercolor = :gold1, markersize = 5);
pl1 = plot!(P_vec2, T2nd_vec2, markershape=:circle, label="Pure shear",lw = 1.5, linecolor =:blue, framestyle = :box, markercolor = :blue, markersize = 5);
pl1 = plot!(P_vec3, T2nd_vec3, markershape=:circle, label="Mixed strain",lw = 1.5, linecolor =:red, framestyle = :box, markercolor = :red, markersize = 5);
pl1 = plot!(P_Y,Tau_y,linecolor =:black,lw = 2, label="", legend=:bottomright, tickfontsize = 13, labelfontsize	= 15, legendfontsize = 12)
pl1 = plot!([Pj, Py],[Tauj, 0], linecolor =:black, label = "", linestyle = :dash,lw = 2)
display(pl1)
savefig("Stress_Pressure.svg")

pl = plot(time_vec2, T2nd_vec2, label=L"τ_\mathrm{II}",lw = 1.5, xlabel="Time [yrs]", ylabel="Stress [MPa]", linecolor =:steelblue1, framestyle = :box);
pl = plot!(time_vec2, P_vec2, label=L"P",lw = 2,legend=:topleft, linecolor =:indianred1, tickfontsize = 13, labelfontsize	= 15, legendfontsize = 12);
display(pl)
savefig("Time_evolution2.svg")

pl = plot(time_vec3, T2nd_vec3, label=L"τ_\mathrm{II}",lw = 1.5, xlabel="Time [yrs]", ylabel="Stress [MPa]", linecolor =:steelblue1, framestyle = :box);
pl = plot!(time_vec3, P_vec3, label=L"P",lw = 2,legend=:topleft, linecolor =:indianred1, tickfontsize = 13, labelfontsize	= 15, legendfontsize = 12);
display(pl)
savefig("Time_evolution3.svg")
=#