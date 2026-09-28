module SEI

using ModelingToolkit
using BatteryToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical

"""
    NoSEI(; name, p::SideReactionParameters, s::SolidParticleParameters, g)

Create a zero SEI growth model (reaction disabled).

Returns a ModelingToolkit system where SEI film thickness remains zero and provides no ohmic resistance.
Use this when SEI growth effects are negligible or you want to exclude them from the simulation.

# Arguments
- `name`: System name for ModelingToolkit (required)
- `p::SideReactionParameters`: SEI reaction parameters (unused in this model)
- `s::SolidParticleParameters`: Electrode solid particle parameters
- `g`: FVM geometry object

# Output Variables
- `L_sei`: SEI film thickness (always 0)
- `j_sei`: SEI current density (always 0)
- `ϕf`: Film potential (always 0)
"""
function NoSEI(; name, params::BatteryParameters, g, Domain)
    
    domain = lowercase(Domain)
    Dmn = split(Domain, " ")[1]  # Extract 'Positive' or 'Negative' for parameter access
    
    @parameters begin
        t
    end

    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    N = g.el.Nx[1]

    @named J = RealInput()
    @named T = RealInput()
    @named Δϕₛ = RealInputArray(nin=N)
    @named aₖ = RealInput(guess=3*params["$(Domain) active material volume fraction"]/params["$(Dmn) particle radius [m]"])

    # Time derivative operator
    Dt = Differential(t)
    
    @variables begin
        # SEI concentration
        (c_sei(t))[1:N] = fill(0,N)
        (j_sei(t))[1:N]
        (ϕf(t))[1:N]
        (L_sei(t))[1:N] = fill(0,N)

        c_sei_x(t)
        L_sei_x(t)
        j_sei_x(t)
        ϕf_x(t)
    end

    η_sei = [Δϕₛ.u[i] - params["SEI open-circuit potential [V]"] - ϕf[i] for i in 1:N]
    
    eqns = [
        # Scott Marquis thesis (eq. 5.92)
        # Exchange current density
        [j_sei[i] ~ 0 for i in 1:N]...,

        [Dt(c_sei[i]) ~ 0 for i in 1:N]...,
        [L_sei[i] ~ c_sei[i]*params["SEI partial molar volume [m3.mol-1]"]/aₖ.u for i in 1:N]...,

        [ϕf[i] ~ J.u*L_sei[i]*params["SEI resistivity [Ohm.m]"] for i in 1:N]...,
        L_sei_x ~ sum(L_sei)/N,
        c_sei_x ~ sum(c_sei)/N,
        j_sei_x ~ sum(j_sei)/N,
        ϕf_x ~ sum(ϕf)/N,
        Q_sei ~ 0,

    ]

    System(eqns,t; name=name,systems=[J, T, Δϕₛ, aₖ])
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
function ReactionLimitedSEI(; name, params::BatteryParameters, g, Domain)
    
    domain = lowercase(Domain)
    Dmn = split(Domain, " ")[1]  # Extract 'Positive' or 'Negative' for parameter access

    @parameters begin
        t
    end

    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    N = g.el.Nx[1]

    @named J = RealInput()
    @named T = RealInput()
    @named Δϕₛ = RealInputArray(nin=N)
    @named aₖ = RealInput(guess=3*params["$(Domain) active material volume fraction"]/params["$(Dmn) particle radius [m]"])

    # Time derivative operator
    Dt = Differential(t)
    
    @variables begin
        # SEI concentration
        (c_sei(t))[1:N] = fill(0,N)#scale
        (j_sei(t))[1:N]
        (ϕf(t))[1:N]
        (L_sei(t))[1:N]

        c_sei_x(t)
        L_sei_x(t)
        j_sei_x(t)
        ϕf_x(t)
        Q_sei(t)
    end

    η_sei = [Δϕₛ.u[i] - params["SEI open-circuit potential [V]"] - ϕf[i] for i in 1:N]
    c_sei₀ = params["Initial SEI thickness [m]"]/params["SEI partial molar volume [m3.mol-1]"]*3*params["$(Domain) active material volume fraction"]/params["$(Dmn) particle radius [m]"]

    eqns = [
        # Scott Marquis thesis (eq. 5.92)
        # Exchange current density
        [j_sei[i] ~ -params["SEI reaction exchange current density [A.m-2]"]*exp(-params["SEI transfer coefficient"]*F/R/T.u*η_sei[i]) for i in 1:N]...,

        [Dt(c_sei[i]) ~ -aₖ.u*j_sei[i]/(F*params["Ratio of lithium moles to SEI moles"]) for i in 1:N]...,
        [L_sei[i] ~ c_sei[i]*params["SEI partial molar volume [m3.mol-1]"]/aₖ.u for i in 1:N]...,

        [ϕf[i] ~ J.u*L_sei[i]*params["SEI resistivity [Ohm.m]"] for i in 1:N]...,
        L_sei_x ~ sum(L_sei)/N,
        c_sei_x ~ sum(c_sei)/N,
        j_sei_x ~ sum(j_sei)/N,
        ϕf_x ~ sum(ϕf)/N,
        Q_sei ~ (c_sei_x-c_sei₀)*params["SEI partial molar volume [m3.mol-1]"]*params["Ratio of lithium moles to SEI moles"]*F/3600,
    ]

    System(eqns,t; name=name,systems=[J, T, Δϕₛ, aₖ])
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
function SolventDiffusionLimitedSEI(; name, p::BatteryToolkit.SideReactionParameters, s::BatteryToolkit.SolidParticleParameters,V, g)
    
    @parameters begin
        t
        D_sol = 2.5e-22
        c_sol = 2636.0, [tunable=false]
        E_sei = 38000.0, [tunable=false]
        T_ref = 298.15
    end

    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    N = g.el.Nx[1]

    @named J = RealInput()
    @named T = RealInput(guess=298.15)
    @named Δϕₛ = RealInputArray(nin=N)
    @named aₖ = RealInput(guess=s.aₖ)

    # Time derivative operator
    Dt = Differential(t)

    c_sei₀ = p.Lf₀/p.V̄*s.aₖ

    @variables begin
        # SEI concentration
        (c_sei(t))[1:N] = fill(c_sei₀,N)
        (j_sei(t))[1:N]
        (ϕf(t))[1:N]
        (L_sei(t))[1:N]

        c_sei_x(t)
        L_sei_x(t)
        j_sei_x(t)
        ϕf_x(t)
        Q_loss(t)
    end

    # All SEI growth mechanisms assumed to have Arrhenius dependence
    arrhenius = exp(
        E_sei / R * (1 / T_ref - 1 / T.u)
    )
    
    eqns = [
        # Scott Marquis thesis (eq. 5.92)
        # Exchange current density
        [j_sei[i] ~ -D_sol*c_sol*F/L_sei[i]*arrhenius for i in 1:N]...,

        [Dt(c_sei[i]) ~ -aₖ.u*j_sei[i]/(F*p.z) for i in 1:N]...,
        [L_sei[i] ~ c_sei[i]*p.V̄/aₖ.u for i in 1:N]...,

        [ϕf[i] ~ -J.u*L_sei[i]*p.R for i in 1:N]...,
        L_sei_x ~ sum(L_sei)/N,
        c_sei_x ~ sum(c_sei)/N,
        j_sei_x ~ sum(j_sei)/N,
        ϕf_x ~ sum(ϕf)/N,
        Q_loss ~ (c_sei_x-c_sei₀)*V*p.z*F/3600,

    ]

    System(eqns,t; name=name,systems=[J, T, Δϕₛ, aₖ])
end

function ECReactionLimitedSEI(; name, p::BatteryToolkit.SideReactionParameters, s::BatteryToolkit.SolidParticleParameters,V, g)
    
    @parameters begin
        t
        k_sei = 2.76e-18
        D_ec = 1.75e-19
        α = 0.5
    end

    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    N = g.el.Nx[1]

    @named J = RealInput()
    @named T = RealInput()
    @named Δϕₛ = RealInputArray(nin=N)
    @named aₖ = RealInput(guess=s.aₖ)

    # Time derivative operator
    Dt = Differential(t)

    c_sei₀ = p.Lf₀/p.V̄*s.aₖ
    c_ec_0 = 4541.0
    # D_ec = 1.75e-19
    # k_sei = 2.76e-18

    @variables begin
        # EC concentration
        (c_ec(t))[1:N], [guess=ones(N)*c_ec_0]
        # SEI concentration
        (c_sei(t))[1:N] = fill(c_sei₀,N)
        (j_sei(t))[1:N]
        (aj_sei(t))[1:N]
        (ϕf(t))[1:N], [guess=zeros(N)]
        (L_sei(t))[1:N], [guess=ones(N)*p.Lf₀]

        c_sei_x(t)
        c_ec_x(t)
        L_sei_x(t), [guess=p.Lf₀]
        j_sei_x(t), [guess=0]
        aj_sei_x(t), [guess=0]
        ϕf_x(t), [guess=0]
        Q_loss(t)
    end

    # All SEI growth mechanisms assumed to have Arrhenius dependence
    arrhenius = exp(
        p.E_sei / R * (1 / p.T_ref - 1 / T.u)
    )

    η_sei = [Δϕₛ.u[i] - p.U + ϕf[i] for i in 1:N]

    k_exp = [k_sei*exp(-α*F/R/T.u*η_sei[i]) for i in 1:N]
    L_over_D = [L_sei[i]/D_ec for i in 1:N]

    eqns = [
        # Scott Marquis thesis (eq. 5.92)
        # Exchange current density
        [j_sei[i] ~ -F*c_ec_0*k_exp[i]/(1 + L_over_D[i]*k_exp[i])*arrhenius for i in 1:N]...,
        [aj_sei[i] ~ aₖ.u*j_sei[i] for i in 1:N]...,
        [c_ec[i] ~ c_ec_0/(1 + L_over_D[i]*k_exp[i]) for i in 1:N]...,

        [Dt(c_sei[i]) ~ -aₖ.u*j_sei[i]/(F*p.z) for i in 1:N]...,
        [L_sei[i] ~ c_sei[i]*p.V̄/aₖ.u for i in 1:N]...,

        [ϕf[i] ~ -J.u*L_sei[i]*p.R for i in 1:N]...,
        L_sei_x ~ sum(L_sei)/N,
        c_sei_x ~ sum(c_sei)/N,
        c_ec_x ~ sum(c_ec)/N,
        j_sei_x ~ sum(j_sei)/N,
        aj_sei_x ~ sum(aj_sei)/N,
        ϕf_x ~ sum(ϕf)/N,
        Q_loss ~ (c_sei_x-c_sei₀)*V*p.z*F/3600,

    ]

    System(eqns,t; name=name,systems=[J, T, Δϕₛ, aₖ])
end

end