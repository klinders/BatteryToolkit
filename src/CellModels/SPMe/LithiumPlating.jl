module LithiumPlating

using ModelingToolkit
using BatteryToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical

"""
    NoPlating(; name, p::SideReactionParameters, s::SolidParticleParameters, g)

Create a zero lithium plating model (reaction disabled).

Returns a ModelingToolkit system where lithium plating is disabled. Use this when lithium
plating effects are negligible or you want to exclude them from the simulation.

# Arguments
- `name`: System name for ModelingToolkit (required)
- `p::SideReactionParameters`: Lithium plating reaction parameters (unused)
- `s::SolidParticleParameters`: Electrode solid particle parameters
- `g`: FVM geometry object

# Output Variables
- `c_plating`: Plated lithium concentration (always 0)
- `c_dead`: Dead lithium concentration (always 0)
- `j_stripping`: Lithium stripping current (always 0)
"""
function NoPlating(; name, p::BatteryToolkit.SideReactionParameters, s::BatteryToolkit.SolidParticleParameters, g)
    
    @parameters begin
        t
    end

    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    N = g.el.Nx[1]

    @named J = RealInput()
    @named T = RealInput()
    @named Δϕₛ = RealInputArray(nin=N)
    @named η_sei = RealInputArray(nin=N)
    @named aₖ = RealInput(guess=s.aₖ)
    @named cₑ = RealInputArray(nin=N)

    # Time derivative operator
    Dt = Differential(t)

    k_plating = 1e-9
    j_strip0 = F*k_plating*1000

    @variables begin
        # Plating concentration
        (c_plating(t))[1:N] = fill(0,N)
        (c_dead(t))[1:N] = fill(0,N)
        (j_stripping(t))[1:N], [guess=j_strip0]
        (ϕf(t))[1:N]

        c_plating_x(t)
        c_dead_x(t)
        j_stripping_x(t)
        ϕf_x(t)

        Q_loss(t)
    end

    eqns = [
        # Scott Marquis thesis (eq. 5.92)
        # Exchange current density
        [j_stripping[i] ~ 0 for i in 1:N]...,
        
        # Irreversable
        [Dt(c_dead[i]) ~ 0 for i in 1:N]...,
        [Dt(c_plating[i]) ~ 0 for i in 1:N]...,
        c_plating_x ~ sum(c_plating)/N,
        c_dead_x ~ sum(c_dead)/N,
        j_stripping_x ~ sum(j_stripping)/N,
        Q_loss ~ 0
    ]

    System(eqns,t; name=name,systems=[J, T, Δϕₛ,η_sei, aₖ, cₑ])
end

"""
    IrreversiblePlating(; name, p::SideReactionParameters, s::SolidParticleParameters, g)

Create an irreversible lithium plating model.

Models lithium plating where all plated lithium becomes "dead" (electrochemically inactive).
Once plated, lithium cannot be stripped (reversed). Use for conservative simulations where
plating is the primary failure mechanism.

# Arguments
- `name`: System name for ModelingToolkit (required)
- `p::SideReactionParameters`: Lithium plating reaction parameters
- `s::SolidParticleParameters`: Electrode solid particle parameters
- `g`: FVM geometry object

# Key Physics
- Plating occurs when surface potential drops too low (overpotential driven)
- All plated lithium irreversibly becomes dead lithium
- Dead lithium consumes cyclable lithium inventory
- Stripping current is zero (irreversible assumption)

# Output Variables
- `c_plating`: Active plated lithium (accumulates until fully dead)
- `c_dead`: Dead lithium (irreversibly lost from battery)
- `j_stripping`: Always zero (no reversibility)
"""
function IrreversiblePlating(; name, p::BatteryToolkit.SideReactionParameters, s::BatteryToolkit.SolidParticleParameters, V, g)
    
    @parameters begin
        t
    end

    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    N = g.el.Nx[1]

    @named J = RealInput()
    @named T = RealInput()
    @named Δϕₛ = RealInputArray(nin=N)
    @named η_sei = RealInputArray(nin=N)
    @named aₖ = RealInput(guess=s.aₖ)
    @named cₑ = RealInputArray(nin=N)
    @named L_sei = RealInputArray(nin=N)

    # Time derivative operator
    Dt = Differential(t)

    scale = 1000 # c_typical
    α_plating = p.α_plating # Li plating transfer coefficient
    α_stripping = 1 - α_plating
    k_plating = p.k_plating
    V̄ = p.V_plating # Partial molar volume of lithium [m3.mol-1]

    j_strip0 = F*k_plating*1000

    @variables begin
        # Plating concentration
        (c_plating(t))[1:N] = zeros(N)
        (c_dead(t))[1:N] = zeros(N)
        (L_plating(t))[1:N]
        (L_dead(t))[1:N]
        (j_stripping(t))[1:N], [guess=fill(0,N)]
        (aj_stripping(t))[1:N]
        (ϕf(t))[1:N]
        (η_plating(t))[1:N]
        (η_stripping(t))[1:N]
        (j0_plating(t))[1:N]
        (j0_stripping(t))[1:N]
        
        c_plating_x(t)
        c_dead_x(t)
        L_plating_x(t)
        L_dead_x(t)
        j_stripping_x(t)
        aj_stripping_x(t)
        ϕf_x(t)
        η_plating_x(t)
        η_stripping_x(t)
        
        Q_plating(t)
        Q_dead(t)
        Q_loss(t)
    end
    
    eqns = [
       [j0_plating[i] ~ F*k_plating*cₑ.u[i] for i in 1:N]...,
        [η_stripping[i] ~ Δϕₛ.u[i] + η_sei.u[i] for i in 1:N]...,
        [η_plating[i] ~ -η_stripping[i] for i in 1:N]...,

        # Scott Marquis thesis (eq. 5.92)
        # Exchange current density
        [j_stripping[i] ~ -j0_plating[i]*exp(α_plating*F/R/T.u*η_plating[i]) 
        for i in 1:N]...,
        [aj_stripping[i] ~ aₖ.u*j_stripping[i] for i in 1:N]...,
        
        # Partially reversible
        [Dt(c_dead[i]) ~  -aₖ.u*j_stripping[i]/F  for i in 1:N]...,
        [Dt(c_plating[i]) ~ 0 for i in 1:N]...,
        c_plating_x ~ sum(c_plating)/N,
        c_dead_x ~ sum(c_dead)/N,
        j_stripping_x ~ sum(j_stripping)/N,
        [L_plating[i] ~ c_plating[i]*V̄/aₖ.u for i in 1:N]...,
        [L_dead[i] ~ c_dead[i]*V̄/aₖ.u for i in 1:N]...,
        L_plating_x ~ sum(L_plating)/N,
        L_dead_x ~ sum(L_dead)/N,
        η_plating_x ~ sum(η_plating)/N,
        η_stripping_x ~ sum(η_stripping)/N,
        aj_stripping_x ~ sum(aj_stripping)/N,
        
        Q_plating ~ c_plating_x*V*F/3600,
        Q_dead ~ c_dead_x*V*F/3600,
        Q_loss ~ Q_plating + Q_dead
    ]

    System(eqns,t; name=name,systems=[J, T, Δϕₛ,η_sei, aₖ, cₑ, L_sei])
end

"""
    PartiallyReversiblePlating(; name, p::SideReactionParameters, s::SolidParticleParameters, g)

Create a partially reversible lithium plating model.

Models lithium plating where a fraction of plated lithium can be stripped (reversed) during
charging, but some becomes dead due to structural damage or electrolyte reactions.
Use for realistic simulations where partial reversibility and cycling damage occur.

# Arguments
- `name`: System name for ModelingToolkit (required)
- `p::SideReactionParameters`: Lithium plating reaction and reversibility parameters
- `s::SolidParticleParameters`: Electrode solid particle parameters
- `g`: FVM geometry object

# Key Physics
- Plating occurs during overdischarge (high negative overpotential)
- Stripping (reversal) occurs during charging with kinetic limitations
- Some plated lithium irreversibly becomes dead due to cycling stress
- Dead lithium fraction accumulates with plating/stripping cycles

# Output Variables
- `c_plating`: Active plated lithium (can grow or shrink with cycling)
- `c_dead`: Dead lithium (irreversibly lost, increases over time)
- `j_stripping`: Stripping current (depends on plating thickness and potential)

# Cycling Behavior
Repeated plating/stripping cycles increase dead lithium fraction, modeling accelerated
degradation under abuse conditions.
"""
function PartiallyReversiblePlating(; name, p::BatteryToolkit.SideReactionParameters, s::BatteryToolkit.SolidParticleParameters, V, g)
    
    @parameters begin
        t
        γ₀ = 5e-6
        α_plating = p.α_plating # Li plating transfer coefficient
        k_plating = p.k_plating
    end

    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    N = g.el.Nx[1]

    @named J = RealInput()
    @named T = RealInput()
    @named Δϕₛ = RealInputArray(nin=N)
    @named η_sei = RealInputArray(nin=N)
    @named aₖ = RealInput(guess=s.aₖ)
    @named cₑ = RealInputArray(nin=N)
    @named L_sei = RealInputArray(nin=N)

    # Time derivative operator
    Dt = Differential(t)

    scale = 1000 # c_typical
    α_stripping = 1 - α_plating
    V̄ = 1.3e-05 # Partial molar volume of lithium [m3.mol-1]

    j_strip0 = F*k_plating*1000
    
    @variables begin
        # Plating concentration
        (c_plating(t))[1:N] = zeros(N)
        (c_dead(t))[1:N] = zeros(N)
        (L_plating(t))[1:N]
        (L_dead(t))[1:N]
        (j_stripping(t))[1:N], [guess=fill(0,N)]
        (aj_stripping(t))[1:N]
        (ϕf(t))[1:N]
        (η_plating(t))[1:N]
        (η_stripping(t))[1:N]
        (j0_plating(t))[1:N]
        (j0_stripping(t))[1:N]
        
        c_plating_x(t)
        c_dead_x(t)
        L_plating_x(t)
        L_dead_x(t)
        j_stripping_x(t)
        aj_stripping_x(t)
        ϕf_x(t)
        η_plating_x(t)
        η_stripping_x(t)
        
        Q_plating(t)
        Q_dead(t)
        Q_loss(t)
    end

    L_sei_0 = 5e-9
    
    coupling = [c_plating[i]*γ₀*L_sei_0/L_sei.u[i] for i in 1:N]

    eqns = [
        [j0_stripping[i] ~ F*k_plating*c_plating[i] for i in 1:N]...,
        [j0_plating[i] ~ F*k_plating*cₑ.u[i] for i in 1:N]...,
        [η_stripping[i] ~ Δϕₛ.u[i] + η_sei.u[i] for i in 1:N]...,
        [η_plating[i] ~ -η_stripping[i] for i in 1:N]...,

        # Scott Marquis thesis (eq. 5.92)
        # Exchange current density
        [j_stripping[i] ~ 
                    j0_stripping[i]*exp(α_stripping*F/R/T.u*η_stripping[i]) -
                    j0_plating[i]*exp(α_plating*F/R/T.u*η_plating[i]) 
        for i in 1:N]...,
        [aj_stripping[i] ~ aₖ.u*j_stripping[i] for i in 1:N]...,
        
        # Partially reversible
        [Dt(c_dead[i]) ~  coupling[i] for i in 1:N]...,
        [Dt(c_plating[i]) ~ -aₖ.u*j_stripping[i]/F - coupling[i] for i in 1:N]...,
        c_plating_x ~ sum(c_plating)/N,
        c_dead_x ~ sum(c_dead)/N,
        j_stripping_x ~ sum(j_stripping)/N,
        [L_plating[i] ~ c_plating[i]*V̄/aₖ.u for i in 1:N]...,
        [L_dead[i] ~ c_dead[i]*V̄/aₖ.u for i in 1:N]...,
        L_plating_x ~ sum(L_plating)/N,
        L_dead_x ~ sum(L_dead)/N,
        η_plating_x ~ sum(η_plating)/N,
        η_stripping_x ~ sum(η_stripping)/N,
        aj_stripping_x ~ sum(aj_stripping)/N,
        
        Q_plating ~ c_plating_x*V*F/3600,
        Q_dead ~ c_dead_x*V*F/3600,
        Q_loss ~ Q_plating + Q_dead
    ]

    System(eqns,t; name=name,systems=[J, T, Δϕₛ,η_sei, aₖ, cₑ, L_sei])
end

end