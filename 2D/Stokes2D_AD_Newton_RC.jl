# 2D Stokes code (V-P formulation) with Newton solver for nonlinearities and AD
using ForwardDiff, GLMakie, SparseArrays, LinearAlgebra
using SparsityTracing, SparseDiffTools
using Test

# Compute strainrate in an allocation free manner
function compute_strainrate!(Exx::Array{_T,2},Ezz::Array{_T,2},Exz::Array{_T,2}, Vx::Array{_T,2},Vz::Array{_T,2}, Δx, Δz) where _T
    diff!(Exx, Vx, Δx, 1)       #Exx .= diff(Vx,dims=1)./Δx;    
    diff!(Ezz, Vz, Δz, 2)       #Ezz .= diff(Vz,dims=2)./Δz;

    for I in eachindex(Exz); Exz[I] = zero(_T); end
    #Exz[2:end-1,2:end-1] .= 0.5.*(diff(Vx[2:end-1,:],dims=2)./Δz + diff(Vz[:,2:end-1],dims=1)./Δx);
    diff!(@view(Exz[2:end-1,2:end-1]), (@view( Vx[2:end-1,:      ]), @view( Vz[:      ,2:end-1])), (2Δz, 2Δx), (2, 1)) 

    return nothing
end

"""
    diff!(diff_a::AbstractArray{T,D}, a::AbstractArray{T,D}, Δ::T, dims::Int=1)

Non-allocating diff function, for a constant grid spacing Δ, along the specified dimension (dims). 
Essentially does `diff_a .= diff(a, dims=dims)./Δ`, but without allocation.
"""
@inline function diff!(diff_a::AbstractArray{T,D}, a::AbstractArray{T,D}, Δ::Number, dims::Int=1) where {T,D}
    one = CartesianIndex(ntuple(i->(i==dims ? 1 : 0), D))
    for I in CartesianIndices(diff_a)
        diff_a[I] = (a[I + one] - a[I]) * inv(Δ)
    end
    nothing
end

# Helper: builds a CartesianIndex with a 1 in position `dim` and 0s elsewhere.
# @generated emits a pure if-else chain with integer literals — no ntuple closure at runtime.
@generated function _diff_one(dim::Int, ::Val{D}) where D
    checks = [:(dim == $d && return CartesianIndex($(ntuple(i -> i==d ? 1 : 0, D)...))) for d in 1:D]
    quote $(checks...) ; CartesianIndex($(ntuple(_ -> 0, D)...)) end
end

"""
    diff!(diff_a::AbstractArray{T,D}, a_tuple::NTuple{N,<:AbstractArray{T,D}}, Δ_tuple::NTuple{N,<:Number}, dims_tuple::NTuple{N,Int})

Non-allocating diff function, for a constant grid spacing Δ, along the specified dimension (dims). 
Essentially does `diff_a .= diff(a[1], dims=dims[1])./Δ[1] + diff(a[2], dims=dims[2])./Δ[2] + ...`, but without allocation.
`diff_a` is zeroed at the start.
"""
@inline @generated function diff!(diff_a::AbstractArray{T,D}, a_tuple::Tuple{Vararg{AbstractArray{T},N}},
               Δ_tuple::NTuple{N,<:Number}, dims_tuple::NTuple{N,Int}) where {T,D,N}
    body = Expr(:block)
    push!(body.args, quote
        @inbounds for I in eachindex(diff_a)
            diff_a[I] = zero(T)
        end
    end)
    for k in 1:N
        push!(body.args, quote
            let A = a_tuple[$k], invΔ = inv(Δ_tuple[$k])
                one = _diff_one(dims_tuple[$k], Val($D))
                @inbounds for I in CartesianIndices(diff_a)
                    diff_a[I] += (A[I + one] - A[I]) * invΔ
                end
            end
        end)
    end
    push!(body.args, :(return nothing))
    return body
end

function power_law_upper_lower_bounded(εII; A=1.0, n=4.5, η_min=1e-3, η_max=1e3)
    # power law with upper and lower bounds
    η_pl = A^(-1/n) * εII^((1-n)/n) 
    η_eff = η_min + inv(inv(η_pl) + inv(η_max))
    
    τII = 2*η_eff*εII

    return τII
end

# This routine calls local routines to update deviatoric stress components
function update_stress(εij, P, phase)
    εII = sqrt(0.5*(εij[1]^2 + εij[2]^2) + εij[3].^2)     # 2nd invariant strainrate tensor

    # to be replaced with RC
    η   = eta_eff_c(εII, phase)
    τij = 2η.*εij
    return τij, η
end

# the rheology routine
function eta_eff_c(εII, phase)
    if phase == 1
        τII = power_law_upper_lower_bounded(εII, n=1)
    else
        τII = power_law_upper_lower_bounded(εII, A=0.001, n=1)
    end

    ηc = τII./(2*εII .+ eps()) # effective viscosity at the center, eps() is added to avoid division by zero
    return ηc
end

function eta_eff_c(εII::AbstractArray, phase_c::AbstractArray)
    ηc = zero(εII)
    for I in eachindex(εII)
        ηc[I] = eta_eff_c(εII[I], phase_c[I])
    end
    return ηc
end

function vertex2center!(center::Array{_T,2}, vertex::Array{_T,2}) where _T
    for I in CartesianIndices(center)
       center[I] =  0.25*(vertex[I[1]+1,I[2]+1] + vertex[I[1]  ,I[2]+1] + vertex[I[1]+1,I[2]  ] + vertex[I[1]  ,I[2]  ])
    end
    return nothing
end

# interpolates from center->vertex; can only be done for vertex[2:end-1,2:end-1]
function center2vertex!(vertex::Array{_T,2}, center::Array{_T,2}) where _T
    nx,nz = size(center)
    for I in CartesianIndices((1:nx-1,1:nz-1))
       vertex[I[1]+1,I[2]+1] =  0.25*(center[I[1]+1,I[2]+1] + center[I[1]  ,I[2]+1] + center[I[1]+1,I[2]  ] + center[I[1]  ,I[2]  ])
    end
    return nothing
end

function second_invariant!(E2nd, Exx, Ezz, Exz2)
    for I in eachindex(E2nd)
       E2nd[I] =  sqrt(0.5*(Exx[I]^2 + Ezz[I]^2) + Exz2[I])
    end
    return nothing
end

# Update stress and pressure
function compute_dev_stress_pressure!(P,Txx,Tzz,Txz,Exx,Ezz,Exz, phase_c, ηc) 
    # Warning: when we deal with VE materials, this must be updated!
    
    # square of the second invariant of the strainrate tensor @ center
    Exz2_c = zero(Exx)
    Exz_c  = zero(Exx)
    Txz_c  = zero(Exx)
    E2nd   = zero(Exx)    
    vertex2center!(Exz2_c, Exz.*Exz)            # interpolate Exz^2 from vertex -> center
    vertex2center!(Exz_c,  Exz)                 # interpolate Exz from vertex -> center
    second_invariant!(E2nd, Exx, Ezz, Exz2_c)   # 2nd invariant

    # update stress components @ center points
    for I in eachindex(Exx) 
        τij, η   = update_stress((Exx[I],Ezz[I], Exz_c[I]), P[I], phase_c[I])  
        Txx[I]   = τij[1]
        Tzz[I]   = τij[2]
        Txz_c[I] = τij[3]
        ηc[I]    = η
    end
   
    # interpolate η to vertex appears more stable than interpolating τxz
    ηv = zero(Exz)
    center2vertex!(ηv, ηc) 
    Txz .= 2ηv.*Exz

    return nothing
end

function compute_residual!(rMass, rVx, rVz, Vx, Vz, P, phase_c, ρc, g, x, z, Δx, Δz, γ, BC)
    nx, nz = size(P)

    ηc  = zeros(eltype(P), nx, nz)
    Exx = zeros(eltype(P), nx, nz)
    Ezz = zeros(eltype(P), nx, nz)
    Exz = zeros(eltype(P), nx+1, nz+1)
    Txx = zeros(eltype(P), nx, nz)
    Tzz = zeros(eltype(P), nx, nz)
    Txz = zeros(eltype(P), nx+1, nz+1)
    
    compute_strainrate!(Exx,Ezz,Exz, Vx,Vz, Δx,Δz)
    compute_dev_stress_pressure!(P, Txx,Tzz,Txz,Exx,Ezz,Exz, phase_c, ηc)

    # Mass conservation residual
    for I in eachindex(rMass)
        rMass[I] = Exx[I] + Ezz[I]  + 1/γ*P[I]
    end
  
    # x-momentum conservation equation
    #rVx[2:end-1,:] .= -diff(P,dims=1)/Δx + diff(Txx,dims=1)/Δx + diff(Txz[2:end-1,:],dims=2)/Δz
    diff!( @view(rVx[2:end-1,:]), (P, Txx,  @view(Txz[2:end-1,:])), (-Δx, Δx, Δz), (1, 1, 2))  

    # z-momentum conservation equation
    #rVz[:,2:end-1] .= -diff(P,dims=2)/Δz + diff(Tzz,dims=2)/Δz + diff(Txz[:,2:end-1],dims=1)/Δx + ρ_vz*g 
    diff!( @view(rVz[:,2:end-1] ), (P, Tzz,  @view(Txz[:,2:end-1])), (-Δz, Δz, Δx), (2, 2, 1))  
    for I in CartesianIndices((1:nx,2:nz-1))
        rVz[I] += 0.5*(ρc[I[1],I[2]] + ρc[I[1],I[2]-1])*g
    end

    # set BC's (zero velocities)
    for iz = 1:nz
        rVx[1,   iz] = Vx[1   , iz] .- x[1  ]*BC.εxx
        rVx[nx+1,iz] = Vx[nx+1, iz] .- x[end]*BC.εxx
    end

    for ix = 1:nz
        rVz[ix,1]    = Vz[ix   ,1   ] .- z[1  ]*BC.εzz
        rVz[ix,nz+1] = Vz[ix   ,nz+1] .- z[end]*BC.εzz
    end

    return nothing
end

# residual routine but in vector format (for use with ForwardDiff)
function residual_vec!(R, X, phase_c, ρc, g, x, z, Δx, Δz, γ, BC)
    nx, nz  = size(phase_c)
    Vx,Vz,P = sol_to_VxVzP(X, nx, nz)

    rVx = zeros(eltype(X), nx+1, nz)
    rVz = zeros(eltype(X), nx, nz+1)
    rP  = zeros(eltype(X), nx, nz)
    compute_residual!(rP, rVx, rVz, Vx, Vz, P, phase_c, ρc, g, x, z, Δx, Δz, γ, BC)
    
    R .= vcat(rVx[:], rVz[:], rP[:])
    return nothing
end

function sol_to_VxVzP(x, nx, nz)
    Vx = reshape(x[1:(nx+1)*nz], nx+1, nz)
    Vz = reshape(x[(nx+1)*nz+1:(nx+1)*nz+(nz+1)*nz], nx, nz+1)
    P  = reshape(x[(nx+1)*nz+(nz+1)*nz+1:end], nx, nz)
    return Vx, Vz, P
end

function set_initial_anomaly!(ρc, xc, zc, inside_value = 2.0, radius=0.1, cen=(0.5,0.5))
    for i in eachindex(xc), j in eachindex(zc)
        if (xc[i]-cen[1])^2 + (zc[j]-cen[2])^2 < radius^2
            ρc[i,j] = inside_value
        end
    end
    return nothing
end


function linesearch(res!, X, ΔX; α=[1e-5,1e-3,1e-2,0.1,0.4,0.9,1.0])
    R, err   = zero(X), zero(α)
    for i in eachindex(α)
        res!(R, X .+ α[i] .* ΔX)
        err[i] = norm(R)
    end
    return α[argmin(err)]
end

"""
    Js, coloring = stokes_solver_sparsity(ρc, phase_c, nx, nz, g, Δx, Δz, γ)
Returns the sparsity pattern of the Jacobian matrix for the 2D Stokes solver. This can be used to optimize the storage and computational efficiency of the Jacobian in the Newton solver.
"""
function stokes_solver_sparsity(ρc, phase_c, nx, nz, g, x, z, γ, BC)
    Δx, Δz  = x[2]-x[1], z[2]-z[1]
 
    # Solution vector
    zc = (z[1:end-1] + z[2:end])/2
    Vx = x[:]* ones(nz)'*(BC.εxx + eps())
    Vz = ones(nx)*z[:]'*(BC.εzz + eps())
    P  = ones(nx)*zc[:]' .* ρc .* g

    # Solution vector
    X       = vcat(Vx[:], Vz[:], P[:])
    X_ad    = SparsityTracing.create_advec(X)
    dR      = similar(X_ad)

    # call residual routine 
    residual_vec!(dR, X_ad, phase_c, ρc, g, x, z, Δx, Δz, γ, BC)

    # Sparsity
    Js      = SparsityTracing.jacobian(dR, length(dR))

    # Coloring
    coloring = SparseDiffTools.matrix_colors(Js)

    return Js, coloring
end


"""
    sol =  stokes_solver_2D(ρc, ηc, ηv, nx, nz, g, Δx, Δz, γ; atol=1e-10, rtol=1e-7, max_it=1500, Js=nothing)

2D nonlinear stokes solver using AD & Newton iterations.
You can optionally provide the sparsity pattern of the Jacobian matrix (Js) along with its coloring to speed up the AD computations. 
If not provided, the Jacobian will be computed as a dense matrix, which can be computationally expensive for large systems.

"""
function stokes_solver_2D(ρc, phase_c, nx, nz, g, x, z, γ, BC; atol=1e-10, rtol=1e-7, max_it=1500, Js=nothing, coloring=nothing)
    Δx, Δz  = x[2]-x[1], z[2]-z[1]

    # Solution vector
    zc = (z[1:end-1] + z[2:end])/2
    Vx = x[:]* ones(nz)'*(BC.εxx + eps())
    Vz = ones(nx)*z[:]'*(BC.εzz + eps())
    P  = ones(nx)*zc[:]' .* ρc .* g

    # Solution vector
    X  = vcat(Vx[:], Vz[:], P[:])
    n  = length(X)
    R  = zeros(n)

    # create an anonymous function to work with ForwardDiff
    res! = (R,X) -> residual_vec!(R, X, phase_c, ρc, g, x, z, Δx, Δz, γ, BC)

    # preallocate Jacobian if possible
    if Js==nothing
        sparsity_provided = false
        J  = zeros(n,n)
        Js = sparse(J)
    else
        sparsity_provided = true
        # build cache for colored Jacobian
        jaccache = ForwardColorJacCache(res!,X;
                              colorvec=coloring,
                              sparsity = Js)
    end

    # inner nonlinear iterations    
    res!(R,X)
    it   = 0
    err  = 1.0
    err0 = norm(R) 
    while err > atol && err>(rtol*err0) && it < max_it 
        if sparsity_provided
            forwarddiff_color_jacobian!(Js, res!, X, jaccache)
        else
            ForwardDiff.jacobian!(J, res!, R, X)
            Js .= sparse(J)
        end
        res!(R,X)
        ΔX   = Js\-R
        α    = linesearch(res!, X, ΔX);
        X    = X + α*ΔX
                
        # compute error 
        err     = norm(R) # requires: using LinearAlgebra 
        it      += 1
        
        println("  Nonlinear iteration $it has error $err, α=$α")
    end
    
    # extract solution
    Vx,Vz,P = sol_to_VxVzP(X, nx, nz)
    Vx_c = (Vx[2:end,:] + Vx[1:end-1,:])/2
    Vz_c = (Vz[:,2:end] + Vz[:,1:end-1])/2

    # postprocessing - compute strainrate
    Exx     = zeros(eltype(P), nx, nz)
    Ezz     = zeros(eltype(P), nx, nz)
    Exz_c   = zeros(eltype(P), nx, nz)
    Exz     = zeros(eltype(P), nx+1, nz+1)
    Eii     = zeros(eltype(P), nx, nz)
    Exz2_c  = zero(Exx)
    compute_strainrate!(Exx,Ezz,Exz, Vx,Vz, Δx,Δz)
    vertex2center!(Exz2_c, Exz.*Exz)            # interpolate Exz^2 from vertex -> center
    second_invariant!(Eii, Exx, Ezz, Exz2_c)    # 2nd invariant


    Txx = zeros(eltype(P), nx, nz)
    Tzz = zeros(eltype(P), nx, nz)
    Tii = zeros(eltype(P), nx, nz)
    ηc  = zeros(nx,nz)
    Txz = zeros(eltype(P), nx+1, nz+1)
    compute_dev_stress_pressure!(P, Txx,Tzz,Txz,Exx,Ezz,Exz, phase_c, ηc)
    Txz2_c  = zero(Txx)
    vertex2center!(Txz2_c, Txz.*Txz)            # interpolate from vertex -> center
    second_invariant!(Tii, Txx, Tzz, Txz2_c)    # 2nd invariant

    return (;Vx, Vz, P, Vx_c, Vz_c, Js, Exx, Ezz, Exz, Eii, Txx,Tzz,Txz,Tii)
end


# Solution grids
nx,nz   = 51, 51
P       = rand(nx  ,nz  )
Vx      = rand(nx+1,nz  )
Vz      = rand(nx  ,nz+1)
Vx[1,:] .= 0.0; Vx[end,:] .= 0.0
Vz[:,1] .= 0.0; Vz[:,end] .= 0.0

# Dimensions
L,H     =   1.0,1.0
Δx,Δz   =   L/nx, H/nz  
x,z     =   0:Δx:L, 0:Δz:H
xc,zc   =   (x[2:end] + x[1:end-1])/2, (z[2:end] + z[1:end-1])/2
εxx     =    1
εzz     =   -1
BC      =   (;εxx,εzz)

# viscosity and density:
phase_c = ones(Int, nx, nz)
ρc      = fill(1.0, nx, nz)
g       = -1.0
ηin     = 1.0
ρin     = 2.0
γ       = 1e8

# Set anomaly of density & viscosity
set_initial_anomaly!(ρc,      xc, zc,  2.0)
set_initial_anomaly!(phase_c, xc, zc,  2)

Js, coloring = stokes_solver_sparsity(ρc, phase_c, nx, nz, g, x, z, γ, BC)

sol = stokes_solver_2D(ρc, phase_c, nx, nz, g, x, z, γ, BC; Js=Js, coloring=coloring)



# Verify solution (for 51x51 case)
@test sum(sol.Vx)   ≈ 1326.0
@test sum(sol.Vz)   ≈ -1326.0000022752758
#@test sum(sol.P)    ≈ -2.4424684852419887e-7
@test sum(sol.Eii)  ≈ 2619.2987615412057
@test sum(sol.Tii)  ≈ 5906.323344013477


# plot
#fig,ax,hm = heatmap(x,z,sol.Exz, axis=(xlabel="x", ylabel="z", title="Exz"))
#fig,ax,hm = heatmap(xc,zc,sol.Exx, axis=(xlabel="x", ylabel="z", title="Exx"))
#fig,ax,hm = heatmap(xc,zc,sol.Tii, axis=(xlabel="x", ylabel="z", title="τII"))
fig,ax,hm = heatmap(xc,zc,sol.Eii, axis=(xlabel="x", ylabel="z", title="εII"))

Colorbar(fig[1, 2], hm)
display(fig)