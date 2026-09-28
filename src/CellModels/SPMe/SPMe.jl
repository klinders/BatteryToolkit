
using ModelingToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical

function abort!(mod,obs,ctx,int)
    ModelingToolkit.terminate!(int)
    t = round(int.t,digits=2)
    @warn "Simulation step terminated at t=$t"
    return (;)
end

include("SolidParticle.jl")
include("Electrolyte.jl")
include("Potentials.jl")
include("SEI.jl")
include("LithiumPlating.jl")
include("ParticleCracking.jl")
include("LAM.jl")
include("CathodeDissolution.jl")

"""
The SPMe model implemented based on [MarquisEtAl2019](@citet) and [BrosaPlanellaWidanage2023](@citet)

**Arguments**
- `name` (optional) defaults to SPMe
- `params ::BatteryParameters` Parameters from the given parameter set
- `Q ::Real` (optional) The capacity of the cell in Ah
- `N ::Dict(:Nₓ=>[::Int, ::Int, ::Int], :Nᵣ=>[::Int, ::Int])` (optional) number of mesh nodes in the particles and electrolyte
- `side_reactions ::Bool` (optional) enable side reactions

"""
function SPMe(; name="SPMe", params::BatteryParameters, Q=0, N=Dict(:Nₓ=>[10,10,10], :Nᵣ=>[10,10]), side_reactions=true)
    @parameters begin
        t # Time variable
    end
    Dt = Differential(t)

    g = build_fvm_geometry(params, N)
    
    if Q==0
        Q = params["Nominal cell capacity [A.h]"]
    end

    Nn = length(g.el.ixₙ)
    Np = length(g.el.ixₚ)

    # Scale the current density to the electrode area
    A = params["Electrode height [m]"]*params["Electrode width [m]"]*params["Number of electrodes connected in parallel to make a cell"]*(Q/params["Nominal cell capacity [A.h]"])
    params["Current collector area [m2]"] = A

    # Electrical ports
    @named p = Pin()
    @named n = Pin()
    @named T = RealInput(guess=298.15)

    # @named Q = RealOutput()
    @named pe = SolidParticle(params, g, d="Positive electrode")
    @named ne = SolidParticle(params, g, d="Negative electrode")
    @named el = Electrolyte(params, g, d=["Negative electrode", "Separator", "Positive electrode"])
    @named sei = SEI.ECReactionLimitedSEI(params, g, d="Negative electrode")
    @named plating = LithiumPlating.PartiallyReversiblePlating(params, g, d="Negative electrode")
    @named cracking_n = ParticleCracking.SwellingAndCracking(params, g, d="Negative electrode")
    @named cracking_p = ParticleCracking.SwellingOnly(params, g, d="Positive electrode")

    @named lam_n = LAM.StressDriven(params. g, d="Negative electrode")
    @named lam_p = LAM.StressDriven(params, g, d="Positive electrode")
    @named cathode_diss = CathodeDissolution.DiffusionCurrent(params, g, d="Positive electrode")

    submodels = [p,n,T,pe,ne,el,sei,plating,cracking_n,cracking_p,lam_n,lam_p,cathode_diss]
    
    cₚ0 = params["Initial concentration in positive electrode [mol.m-3]"]
    cₙ0 = params["Initial concentration in negative electrode [mol.m-3]"]
    cₙ₊ = params["Maximum concentration in negative electrode [mol.m-3]"]
    cₚ₊ = params["Maximum concentration in positive electrode [mol.m-3]"]

    Up_init = params["Positive electrode OCP [V]"](cₚ0/cₚ₊)
    Un_init = params["Negative electrode OCP [V]"](cₙ0/cₙ₊)

    @variables begin
        # Terminal voltage and current
        v(t), [guess=Up_init - Un_init]
        i(t)
        soc(t)

        U₀(t)
        ηᵣ(t)
        ηᵣ̅n(t)
        ηᵣ̅p(t)
        (ηₙ(t))[1:Nn]
        (ηₚ(t))[1:Np]
        Δϕₛ(t), [guess=0]
        (Δϕₙ(t))[1:Nn], [guess=fill(Un_init, Nn)]
        (Δϕₚ(t))[1:Np], [guess=fill(Up_init, Np)]

        # Δϕf(t)
        Δϕₙ_x(t)
        Δϕₙ_x2(t)
        (ϕₙ(t))[1:g.el.Nx[1]]
        (ϕₚ(t))[1:g.el.Nx[3]]
        (jₙ0(t))[1:g.el.Nx[1]]
        (jₚ0(t))[1:g.el.Nx[3]]
        j̄ₙ0(t)
        j̄ₚ0(t)
        ϕ̄ₙ(t)
        ϕ̄ₚ(t)
        Cₙ(t) # Negative electrode capacity in Ah
        Cₚ(t) # Positive electrode capacity in Ah
        Q_loss(t)
        C_cell(t)
        C_pos(t)
        C_neg(t)
        Q_Ah(t) = 0
        Qt_Ah(t) = 0
        n_Li(t)
        LAMₚ(t)
        LAMₙ(t)
        LLI(t)

        #temp
        j_tot_ne(t)
        aj_tot_ne(t)

        Rᵢ(t)
    end
    
    i_app = i/A

    # X-average
    x = g.el.x_centers
    Lₙ = params["Negative electrode thickness [m]"]
    Lₚ = params["Positive electrode thickness [m]"]
    Lₛ = params["Separator thickness [m]"]
    L = Lₙ + Lₛ + Lₚ

    σₙ = params["Negative electrode conductivity [S.m-1]"]
    σₚ = params["Positive electrode conductivity [S.m-1]"]

    # Exchange current densities

    ## Reaction overpotentials ##

    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    
    # Overpotentials from the inverse Butler-Volmer equation
    ηᵣn = [2*R*T.u/F*asinh((i_app/Lₙ/ne.aₖ)/(2*jₙ0[i])) for i in 1:Nn]
    ηᵣp = [2*R*T.u/F*asinh((-i_app/Lₚ/pe.aₖ)/(2*jₚ0[i])) for i in 1:Np]

    # ηᵣ_n = ηᵣ_x(params.n, cₙ, el.cₑ[g.el.ixₙ], ne.T.u, ne.J.u)
    # ηᵣ_p = ηᵣ_x(params.p, cₚ, el.cₑ[g.el.ixₚ], pe.T.u, pe.J.u)

    # X-average of the electrolyte potential
    ϕₛ_n = [i_app*(x[i] - 2*Lₙ)*x[i]/2/σₙ/Lₙ for i in g.el.ixₙ]
    ϕₛ_p = [i_app*(x[i] + (x[i] - L)^2/(2*Lₚ))/σₚ for i in g.el.ixₚ]

    eqns = [
        # Temps
        el.T.u ~ T.u,
        pe.T.u ~ T.u,
        ne.T.u ~ T.u,
        soc ~ ne.z,

        # temp
        aj_tot_ne ~ i_app/Lₙ,
        j_tot_ne ~ aj_tot_ne/ne.aₖ,

        ## Potentials ##
        U₀ ~ pe.U₀ - ne.U₀,
        ηᵣ̅n ~ sum(ηᵣn)/Nn,
        ηᵣ̅p ~ sum(ηᵣp)/Np,
        ηᵣ ~ ηᵣ̅p - ηᵣ̅n,
        Δϕₛ ~ -i_app/3*(Lₚ/σₚ + Lₙ/σₙ),

        # Exchange current densities
        [jₙ0[i] ~ params["Negative electrode exchange current density [A.m-2]"](el.cₑ[g.el.ixₙ[i]], ne.c_surf, params["Maximum concentration in negative electrode [mol.m-3]"], T.u) for i in 1:Nn]...,
        [jₚ0[i] ~ params["Positive electrode exchange current density [A.m-2]"](el.cₑ[g.el.ixₚ[i]], pe.c_surf, params["Maximum concentration in positive electrode [mol.m-3]"], T.u) for i in 1:Np]...,
        j̄ₙ0 ~ sum(jₙ0)/Nn,
        j̄ₚ0 ~ sum(jₚ0)/Np,

        [ϕₙ[i] ~ n.v + ϕₛ_n[i] for i in 1:Nn]...,
        [ϕₚ[i] ~ p.v - ϕₛ_p[i] for i in 1:Np]...,
        [ηₙ[i] ~ ϕₙ[i] - el.ϕₑ[g.el.ixₙ[i]] for i in 1:Nn]...,
        [ηₚ[i] ~ ϕₚ[i] - el.ϕₑ[g.el.ixₚ[i]] for i in 1:Np]...,
        ϕ̄ₙ ~ sum(ϕₙ)/Nn,
        ϕ̄ₚ ~ sum(ϕₚ)/Np,
        v ~ U₀ + ηᵣ + el.ηₑ + el.Δϕₑ + Δϕₛ + sei.ϕf_x,
        Rᵢ ~ el.Δϕₑ + Δϕₛ + sei.ϕf_x, 

        v ~ p.v - n.v,
        p.i ~ i,
        n.i ~-i,

        # Potential differences
        [Δϕₙ[i] ~ ϕₙ[i] - el.ϕₑ[g.el.ixₙ[i]] for i in 1:Nn]...,
        [Δϕₚ[i] ~ ϕₚ[i] - el.ϕₑ[g.el.ixₚ[i]] for i in 1:Np]...,
        Δϕₙ_x ~ ne.U₀ + ϕ̄ₙ + ηᵣ̅n - sei.ϕf_x,

        # Electrolyte current density
        el.i_app.u ~ i_app,
        el.j_n.u ~ i_app/Lₙ,
        el.j_p.u ~ -i_app/Lₚ,
        # el.jₙ0.u ~ jₙ0,
        el.ϕₛn.u ~ ϕ̄ₙ,
        el.Δϕₙ.u ~ ne.U₀ + ηᵣ̅n - sei.ϕf_x,

        ne.J.u ~  i_app/Lₙ/ne.aₖ - sei.j_sei_x - plating.j_stripping_x - cracking_n.j_sei_x, # Current density in the negative electrode
        pe.J.u ~  -i_app/Lₚ/pe.aₖ, # Current density in the positive electrode

        # # Ne sei reaction
        sei.J.u ~ ne.J.u, # Total current density for SEI potential
        sei.T.u ~ T.u,
        sei.aₖ.u ~ ne.aₖ,
        [sei.Δϕₛ.u[i] ~ Δϕₙ_x for i in 1:Nn]...,
        sei.k_sei ~ k_sei,
        sei.D_ec ~ D_ec,
        sei.α ~ α,
        cracking_n.k_sei ~ k_sei,
        cracking_n.D_ec ~ D_ec,
        cracking_n.α ~ α,
        # sei.D_sol ~ D_sol,

        # # Li plating
        plating.J.u ~ i_app/Lₙ/ne.aₖ,
        plating.T.u ~ T.u,
        plating.aₖ.u ~ ne.aₖ,
        [plating.Δϕₛ.u[i] ~ Δϕₙ_x for i in 1:Nn]...,
        plating.η_sei.u ~ sei.ϕf,
        [plating.cₑ.u[i] ~ el.cₑ[i] for i in 1:Nn]...,
        [plating.L_sei.u[i] ~ sei.L_sei[i] for i in 1:Nn]...,

        ## Cracking
        cracking_n.J.u ~ ne.J.u,
        cracking_n.T.u ~ T.u,
        cracking_n.aₖ.u ~ ne.aₖ,
        [cracking_n.Δϕₛ.u[i] ~ Δϕₙ_x for i in 1:Nn]...,
        cracking_n.c_s_r.u ~ ne.c_r,
        cracking_n.c_s_surf.u ~ ne.c_surf,

        cracking_p.J.u ~ -i_app/Lₚ/pe.aₖ,
        cracking_p.T.u ~ T.u,
        cracking_p.aₖ.u ~ pe.aₖ,
        [cracking_p.Δϕₛ.u[i] ~ ϕₚ[i] + el.ϕₑ[g.el.ixₚ[i]] for i in 1:Np]...,
        cracking_p.c_s_r.u ~ pe.c_r,
        cracking_p.c_s_surf.u ~ pe.c_surf,

        ## LAM
        lam_n.σₜ.u ~ cracking_n.σₜ,
        lam_n.σᵣ.u ~ cracking_n.σᵣ,
        lam_n.c_r.u ~ ne.c_r,
        Dt(ne.ϵₛ) ~ lam_n.j_lam,
        
        lam_p.σₜ.u ~ cracking_p.σₜ,
        lam_p.σᵣ.u ~ cracking_p.σᵣ,
        lam_p.c_r.u ~ pe.c_r,
        cathode_diss.Δϕₛ.u ~ Δϕₙ_x,
        cathode_diss.T.u ~ T.u,
        Dt(pe.ϵₛ) ~ lam_p.j_lam + cathode_diss.j_diss,
        
        # Porosity (assumed constant)
        [el.ϵ[i] ~ params["Negative electrode porosity"] - ne.aₖ*(
            sei.L_sei[i] - params["Initial SEI thickness [m]"]
            + plating.L_plating[i] 
            + plating.L_dead[i]
             + cracking_n.L_sei[i]*(cracking_n.r_surf - 1)
            ) for i in g.el.ixₙ]...,
        [el.ϵ[i] ~ params["Separator porosity"] for i in g.el.ixₛ]...,
        [el.ϵ[i] ~ params["Positive electrode porosity"] for i in g.el.ixₚ]...,

        Cₚ ~ pe.ϵₛ*Lₚ * A * params["Maximum concentration in positive electrode [mol.m-3]"] * F / 3600,
        Cₙ ~ ne.ϵₛ*Lₙ * A * params["Maximum concentration in negative electrode [mol.m-3]"] * F / 3600,

        Q_loss ~ sei.Q_loss + plating.Q_loss + cracking_n.Q_sei + lam_n.Q_loss + lam_p.Q_loss,

        C_cell ~ Q - Q_loss,
        C_pos ~ Lₚ*A*params["Maximum concentration in positive electrode [mol.m-3]"]*pe.ϵₛ*F/3600,
        C_neg ~ Lₙ*A*params["Maximum concentration in negative electrode [mol.m-3]"]*ne.ϵₛ*F/3600,
        n_Li ~ sum(el.cₑ) + sum(ne.c) + sum(pe.c),

        LAMₚ ~ (1- C_pos/(Lₚ*A*params["Maximum concentration in positive electrode [mol.m-3]"]*params["Positive electrode active material volume fraction"]*F/3600))*100,
        LAMₙ ~ (1- C_neg/(Lₙ*A*params["Maximum concentration in negative electrode [mol.m-3]"]*params["Negative electrode active material volume fraction"]*F/3600))*100,
        LLI ~ (1 - n_Li/(params["Initial concentration in electrolyte [mol.m-3]"]*g.el.Nₜ + params["Initial concentration in negative electrode [mol.m-3]"]*Nn + params["Initial concentration in positive electrode [mol.m-3]"]*Np))*100,

        Dt(Q_Ah) ~ i/3600,
        Dt(Qt_Ah) ~ abs(i)/3600

    ]

    # Event working
    events = [
        [
            v ~ params["Lower voltage cut-off [V]"],
            v ~ params["Upper voltage cut-off [V]"],
        ]=>(abort!,(;))
    ]

    return System(eqns, t; name=name,systems=submodels, continuous_events=events)
end