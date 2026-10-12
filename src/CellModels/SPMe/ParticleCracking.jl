module ParticleCracking

using ModelingToolkit
using BatteryToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical

"""
    NoCracking(; name, p::SideReactionParameters, s::SolidParticleParameters, g)

Create a model without cracking (reaction disabled).

Returns a ModelingToolkit system where crack thickness and SEI layer remains zero and provides no ohmic resistance.
Use this when cracking effects are negligible or you want to exclude them from the simulation.

# Arguments
- `name`: System name for ModelingToolkit (required)
- `p::SideReactionParameters`: Cracking reaction parameters (unused in this model)
- `s::SolidParticleParameters`: Electrode solid particle parameters
- `g`: FVM geometry object

# Output Variables
- `L_sei`: SEI film thickness (always 0)
- `j_sei`: SEI current density (always 0)
- `ϕf`: Film potential (always 0)
"""
function NoMechanics(; name, p::BatteryToolkit.SideReactionParameters, s::BatteryToolkit.SolidParticleParameters,V ,g)
    
    @parameters begin
        t
    end

    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    N = g.el.Nx[1]

    @named J = RealInput()
    @named T = RealInput()
    @named Δϕₛ = RealInputArray(nin=N)
    @named aₖ = RealInput(guess=s.aₖ)
    @named c_s_r = RealInput()
    @named c_s_surf = RealInput()


    # Time derivative operator
    Dt = Differential(t)
    

    
    @variables begin
        l_cr(t)
        r_surf(t)
        a_cr(t)
        σₜ(t)
        σᵣ(t)
        u_d(t)
        
        # SEI concentration
        (c_sei(t))[1:N] = fill(0,N)
        (j_sei(t))[1:N]
        (ϕf(t))[1:N]
        (L_sei(t))[1:N]

        l_cr_x(t)
        r_surf_x(t)
        a_cr_x(t)
        σₜ_x(t)
        σᵣ_x(t)
        u_d_x(t)

        c_sei_x(t)
        L_sei_x(t)
        j_sei_x(t)
        ϕf_x(t)
        Q_sei(t)
    end

    # dK_SIF = ifelse(σₜ >= 0, σₜ*b_cr* sqrt(pi*l_cr), 0)

    eqns = [

        # Cracking
        l_cr ~ 0,
        r_surf ~ 1,
        a_cr ~ 0,
        σₜ ~ 0,
        σᵣ ~ 0,
        u_d ~ 0,

        # Scott Marquis thesis (eq. 5.92)
        # Exchange current density
        [j_sei[i] ~ 0 for i in 1:N]...,

        [Dt(c_sei[i]) ~ 0 for i in 1:N]...,
        [L_sei[i] ~ 0 for i in 1:N]...,

        [ϕf[i] ~ 0 for i in 1:N]...,

        l_cr_x ~ l_cr,
        r_surf_x ~ r_surf,
        a_cr_x ~ a_cr,
        σₜ_x ~ σₜ,
        σᵣ_x ~ σᵣ,
        u_d_x ~ u_d,

        L_sei_x ~ sum(L_sei)/N,
        c_sei_x ~ sum(c_sei)/N,
        j_sei_x ~ sum(j_sei)/N,
        ϕf_x ~ sum(ϕf)/N,
        Q_sei ~ 0,
    ]

    System(eqns,t; name=name,systems=[J, T, Δϕₛ, aₖ, c_s_r, c_s_surf])
end


"""
    ReactionLimitedSEI(; name, p::SideReactionParameters, s::SolidParticleParameters, g)

Create a reaction-limited SEI growth model.

Models SEI film formation with reaction kinetics controlled by surface overpotential.
The SEI current density follows Butler-Volmer kinetics. Use when SEI growth is fast
(high overpotential) and film diffusion resistance is negligible.

# Arguments
- `name`: System name for ModelingToolkit (required)
- `p::SideReactionParameters`: SEI reaction kinetic parameters
- `s::SolidParticleParameters`: Electrode solid particle parameters
- `g`: FVM geometry object

# Key Parameters Used
- `p.j_sei₀`: Exchange current density (A/m²)
- `p.α`: Transfer coefficient (charge transfer kinetics)
- `p.U`: SEI formation potential (V vs Li/Li⁺)

# Output Variables
- `L_sei`: SEI film thickness (grows over time)
- `j_sei`: SEI current density (determined by kinetics)
- `ϕf`: Film potential (typically small)

# Physical Assumption
Reaction rate dominates over diffusion; film acts as perfect ionic conductor.
"""
function SwellingOnly(; name, p::BatteryToolkit.SideReactionParameters, s::BatteryToolkit.SolidParticleParameters,V ,g)
    
    @parameters begin
        t
    end

    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    N = g.el.Nx[1]


    # Physical constants for stress/displacement calculations (Ai2019)
    stress_geometric_factor = 3.0
    displacement_geometric_factor = 3.0
    crack_roughness_factor = 2.0

    R0 = s.Rₖ
    c₀ = 0

    @named J = RealInput()
    @named T = RealInput()
    @named Δϕₛ = RealInputArray(nin=N)
    @named aₖ = RealInput(guess=s.aₖ)
    @named c_s_r = RealInput()
    @named c_s_surf = RealInput()


    # Time derivative operator
    Dt = Differential(t)
    
    # All SEI growth mechanisms assumed to have Arrhenius dependence
    arrhenius = exp(
        p.E_sei / R * (1 / p.T_ref - 1 / T.u)
    )
    
    @variables begin
        l_cr(t)
        r_surf(t)
        a_cr(t)
        σₜ(t)
        σᵣ(t)
        u_d(t)
        
        # SEI concentration
        (c_sei(t))[1:N] = fill(0,N)
        (j_sei(t))[1:N]
        (ϕf(t))[1:N]
        (L_sei(t))[1:N]

        l_cr_x(t)
        r_surf_x(t)
        a_cr_x(t)
        σₜ_x(t)
        σᵣ_x(t)
        u_d_x(t)

        c_sei_x(t)
        L_sei_x(t)
        j_sei_x(t)
        ϕf_x(t)
        Q_sei(t)
    end

    # dK_SIF = ifelse(σₜ >= 0, σₜ*b_cr* sqrt(pi*l_cr), 0)

    eqns = [

        # Cracking
        l_cr ~ 0,
        r_surf ~ 1,
        a_cr ~ 0,
        σₜ ~ s.Ω*s.E*(c_s_r.u - c_s_surf.u)/stress_geometric_factor/(1.0 - s.ν),
        σᵣ ~ 0,
        u_d ~ s.Ω*R0*(c_s_r.u - c₀)/displacement_geometric_factor,

        # Scott Marquis thesis (eq. 5.92)
        # Exchange current density
        [j_sei[i] ~ 0 for i in 1:N]...,

        [Dt(c_sei[i]) ~ 0 for i in 1:N]...,
        [L_sei[i] ~ 0 for i in 1:N]...,

        [ϕf[i] ~ 0 for i in 1:N]...,

        l_cr_x ~ l_cr,
        r_surf_x ~ r_surf,
        a_cr_x ~ a_cr,
        σₜ_x ~ σₜ,
        σᵣ_x ~ σᵣ,
        u_d_x ~ u_d,

        L_sei_x ~ sum(L_sei)/N,
        c_sei_x ~ sum(c_sei)/N,
        j_sei_x ~ sum(j_sei)/N,
        ϕf_x ~ sum(ϕf)/N,
        Q_sei ~ 0,
    ]

    System(eqns,t; name=name,systems=[J, T, Δϕₛ, aₖ, c_s_r, c_s_surf])
end

"""
    SolventDiffusionLimitedSEI(; name, p::SideReactionParameters, s::SolidParticleParameters, g)

Create a solvent-diffusion-limited SEI growth model.

Models SEI film formation limited by solvent diffusion through the growing film.
The SEI current density decreases as the film thickens due to increasing ionic resistance
and decreasing solvent diffusion. Use when film resistance dominates over reaction kinetics.

# Arguments
- `name`: System name for ModelingToolkit (required)
- `p::SideReactionParameters`: SEI reaction and film transport parameters
- `s::SolidParticleParameters`: Electrode solid particle parameters
- `g`: FVM geometry object

# Key Parameters Used
- `p.D_sol`: Solvent diffusivity in film (m²/s)
- `p.c_sol`: Solvent concentration (mol/m³)
- `p.U`: SEI formation potential (V vs Li/Li⁺)
- `p.R`: Film resistivity (Ω·m)

# Output Variables
- `L_sei`: SEI film thickness (grows over time, asymptotically)
- `j_sei`: SEI current density (decreases as L_sei increases)
- `ϕf`: Film potential drop (increases with thickness)

# Physical Assumption
Film diffusion resistance and potential drop dominate; SEI growth self-limits via thickness.
"""
function SwellingAndCracking(; name, p::BatteryToolkit.SideReactionParameters, s::BatteryToolkit.SolidParticleParameters, V, g)
    
    @parameters begin
        t
        k_sei = 2.76e-18
        D_ec = 1.75e-19
        α = 0.5
        k_cr = 3.9e-20
        b_cr = 1.12
        m_cr = 2.2
    end

    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    N = g.el.Nx[1]


    # Physical constants for stress/displacement calculations (Ai2019)
    stress_geometric_factor = 3.0
    displacement_geometric_factor = 3.0
    crack_roughness_factor = 2.0

    R0 = s.Rₖ
    c₀ = 0

    @named J = RealInput()
    @named T = RealInput()
    @named Δϕₛ = RealInputArray(nin=N)
    @named aₖ = RealInput(guess=s.aₖ)
    @named c_s_r = RealInput()
    @named c_s_surf = RealInput()


    # Time derivative operator
    Dt = Differential(t)

    l_cr_0 = 2e-8
    L₀ = 5e-13
    c_sei₀ = L₀/p.V̄*s.aₖ
    c_ec_0 = 4541.0
    # D_ec = 1.75e-19
    # k_sei = 2.76e-16
    # D_sol = 2.5e-22
    # c_sol = 2636.0
    
    # All SEI growth mechanisms assumed to have Arrhenius dependence
    arrhenius = exp(
        p.E_sei / R * (1 / p.T_ref - 1 / T.u)
    )
    
    @variables begin
        l_cr(t) = l_cr_0
        r_surf(t)
        a_cr(t)
        σₜ(t)
        σᵣ(t)
        u_d(t)
        
        # SEI concentration
        (c_ec(t))[1:N], [guess=fill(c_ec_0, N)]
        (c_sei(t))[1:N] = fill(c_sei₀,N)
        (j_sei(t))[1:N]
        (aj_sei(t))[1:N]
        (ϕf(t))[1:N]
        (L_sei(t))[1:N]

        l_cr_x(t)
        r_surf_x(t)
        a_cr_x(t)
        σₜ_x(t)
        σᵣ_x(t)
        u_d_x(t)

        c_ec_x(t)
        c_sei_x(t)
        L_sei_x(t)
        j_sei_x(t)
        aj_sei_x(t)
        ϕf_x(t)
        Q_sei(t)
    end

    dK_SIF = ifelse(σₜ >= 0, σₜ*b_cr* sqrt(pi*l_cr), 0)

    η_sei = [Δϕₛ.u[i] - p.U + ϕf[i] for i in 1:N]

    k_exp = [k_sei*exp(-α*F/R/T.u*η_sei[i]) for i in 1:N]
    L_over_D = [L_sei[i]/D_ec for i in 1:N]

    eqns = [

        # Cracking
        Dt(l_cr) ~ k_cr*(dK_SIF^m_cr)/3600,
        r_surf ~ 1 + crack_roughness_factor*l_cr*s.ρ_cr*s.w_cr,
        a_cr ~ (r_surf - 1)*aₖ.u,
        σₜ ~ s.Ω*s.E*(c_s_r.u - c_s_surf.u)/stress_geometric_factor/(1.0 - s.ν),
        σᵣ ~ 0,
        u_d ~ s.Ω*R0*(c_s_r.u - c₀)/displacement_geometric_factor,

        # Scott Marquis thesis (eq. 5.92)
        # Exchange current density
        # Scott Marquis thesis (eq. 5.92)
        # Exchange current density
        # [j_sei[i] ~ -D_sol*c_sol*F/L_sei[i]*arrhenius for i in 1:N]...,

        [j_sei[i] ~ -F*c_ec_0*k_exp[i]/(1 + L_over_D[i]*k_exp[i])*arrhenius for i in 1:N]...,
        [aj_sei[i] ~ a_cr*j_sei[i] for i in 1:N]...,
        [c_ec[i] ~ c_ec_0/(1 + L_over_D[i]*k_exp[i]) for i in 1:N]...,

        [Dt(c_sei[i]) ~ -a_cr*j_sei[i]/(F*p.z) for i in 1:N]...,
        [L_sei[i] ~ c_sei[i]*p.V̄/a_cr for i in 1:N]...,

        [ϕf[i] ~ -J.u*L_sei[i]*p.R for i in 1:N]...,

        l_cr_x ~ l_cr,
        r_surf_x ~ r_surf,
        a_cr_x ~ a_cr,
        σₜ_x ~ σₜ,
        σᵣ_x ~ σᵣ,
        u_d_x ~ u_d,

        L_sei_x ~ sum(L_sei)/N,
        # c_ec_x ~ sum(c_ec)/N,
        c_sei_x ~ sum(c_sei)/N,
        j_sei_x ~ sum(j_sei)/N,
        aj_sei_x ~ sum(aj_sei)/N,
        ϕf_x ~ sum(ϕf)/N,
        Q_sei ~ (c_sei_x-c_sei₀)*V*p.z*F/3600,
    ]

    System(eqns,t; name=name,systems=[J, T, Δϕₛ, aₖ, c_s_r, c_s_surf])
end

end