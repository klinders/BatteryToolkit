module CathodeDissolution

using ModelingToolkit
using BatteryToolkit
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical

function NoDiffusion(; name, s::BatteryToolkit.SolidParticleParameters,V, g)
    
    @parameters begin
        t
    end

    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    N = g.el.Nx[1]

    @named T = RealInput(guess=298.15)
    @named Δϕₛ = RealInput()

    # Time derivative operator
    Dt = Differential(t)

    @variables begin
        j_diss(t)
    end

    eqns = [
        # Scott Marquis thesis (eq. 5.92)
        # Exchange current density
        j_diss ~ 0,

    ]

    System(eqns,t; name=name,systems=[T, Δϕₛ])
end

function DiffusionCurrent(; name, s::BatteryToolkit.SolidParticleParameters,V, g)
    
    @parameters begin
        t
        i₀ = 6.24e−04
        Ediss = 4
    end

    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant
    N = g.el.Nx[1]

    @named T = RealInput(guess=298.15)
    @named Δϕₛ = RealInput()

    # Time derivative operator
    Dt = Differential(t)

    @variables begin
        j_diss(t)
    end

    η_diss = Δϕₛ.u - Ediss

    eqns = [
        # Scott Marquis thesis (eq. 5.92)
        # Exchange current density
        j_diss ~ -i₀*exp(F/(R*T.u)*η_diss),

    ]

    System(eqns,t; name=name,systems=[T, Δϕₛ])
end
    
end