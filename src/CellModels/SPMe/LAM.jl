module LAM

using ModelingToolkit
using BatteryToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical


function NoLAM(; name, p::BatteryToolkit.SideReactionParameters, s::BatteryToolkit.SolidParticleParameters, V, g)
    
    @parameters begin
        t
    end
    
    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    N = g.el.Nx[1]
    
    @named σᵣ = RealInput()
    @named σₜ = RealInput()
    @named c_r = RealInput()

    # Time derivative operator
    Dt = Differential(t)
    
    @variables begin
        j_lam(t)
        lli(t) = 0
        Q_loss(t)
    end

    eqns = [
        j_lam ~ 0,
        Dt(lli) ~ 0,
        Q_loss ~ lli * F / 3600
    ]

    System(eqns,t; name=name,systems=[σᵣ, σₜ, c_r])
end

"""
    NoLAM(; name, p::SideReactionParameters, s::SolidParticleParameters, g)

Create a zero LAM model (reaction disabled).

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
function StressDriven(; name, p::BatteryToolkit.SideReactionParameters, s::BatteryToolkit.SolidParticleParameters, V, g)
    
    @parameters begin
        t
        β_LAM = 2.7778e-07 # LAM rate constant
        m_LAM = 2.0 # LAM stress exponent
    end
    
    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    N = g.el.Nx[1]
    
    @named σᵣ = RealInput()
    @named σₜ = RealInput()
    @named c_r = RealInput()


    # Time derivative operator
    Dt = Differential(t)
    
    @variables begin
        j_lam(t)
        lli(t) = 0
        Q_loss(t)
    end
    
    # obtain the rate of loss of active materials (LAM) by stress
    # This is loss of active material model by mechanical effects
    
    # compute hydrostatic stress
    σₕ = (σᵣ.u + 2 * σₜ.u) / 3

    # separate compressive and tensile stresses
    σₕ_c = ifelse(σₕ < 0 , σₕ, 0)
    σₕ_t = ifelse(σₕ > 0 , σₕ, 0)
    
    # assuming that only tensile stress contributes and that the minimum
    # (tensile) hydrostatic stress is zero for full cycles
    σₕ_min = σₕ * 0

    eqns = [
        j_lam ~ -β_LAM*((σₕ_t - σₕ_min) / s.stress_critical)^m_LAM,
        Dt(lli) ~ -V * c_r.u * j_lam,
        Q_loss ~ lli * F / 3600
    ]

    System(eqns,t; name=name,systems=[σᵣ, σₜ, c_r])
end

end