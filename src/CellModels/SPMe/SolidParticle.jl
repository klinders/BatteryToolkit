using ModelingToolkit

"""
    SolidParticle(; name, p::SolidParticleParameters, g)

Create a ModelingToolkit system for solid-state lithium diffusion in a battery electrode.

Models radial diffusion of lithium ions within spherical electrode particles using the finite
volume method. Computes surface concentration, stoichiometry, and open-circuit potential.

# Arguments
- `name`: System name for ModelingToolkit (required)
- `p::SolidParticleParameters`: Electrode material parameters
- `g`: FVM geometry object with node locations and volumes

# Input Ports
- `J`: Surface current density (A/m²)
- `T`: Temperature (K)

# Output Variables
- `c`: Concentration profile across particles (mol/m³)
- `c_surf`: Surface lithium concentration (mol/m³)
- `z`: Stoichiometry = c_surf/c_max
- `U₀`: Open-circuit potential (V)

# Notes
Uses second-order accurate finite volume discretization with ghost nodes for boundary conditions.
"""
function SolidParticle(; name, params::BatteryParameters, g, Domain)
    
    @assert domain in ["positive electrode", "negative electrode"] "Domain must be either 'positive electrode' or 'negative electrode'"
    
    @parameters begin
        t
    end

    R = 8.314 # Universal gas constant
    F = 96485 # Faraday's constant

    @named J = RealInput()
    @named T = RealInput(guess=298.15)

    # Time derivative operator
    Dt = Differential(t)
    domain = lowercase(Domain)
    Dmn = split(Domain, " ")[1]  # Extract 'Positive' or 'Negative' for parameter access

    @variables begin
        # I am adding two ghost nodes for the boundary conditions
        (c(t))[1:g.Nᵣ] = fill(params["Initial concentration in $(domain) [mol.m-3]"],g.Nᵣ)
        (D(t))[1:g.Nᵣ]
        (σ(t))[1:g.Nᵣ]
        D_r(t)
        c_avr(t)
        c_r(t)
        c_surf(t)
        U₀(t)
        z(t)
        ϕ̄ₛ(t)
        ϵₛ(t) = params["$(Domain) active material volume fraction"] # active material volume fraction
        aₖ(t) # specific surface area
    end

    # Discretized equations
    Δr,r,Vᵢ,Aₗ,Aᵣ = g.Δr, g.r, g.Vᵢ, g.Aₗ, g.Aᵣ

    Ω = params["$(Domain) partial molar volume [m3.mol-1]"]
    ν = params["$(Domain) Poisson's ratio"]
    E = params["$(Domain) Young's modulus [Pa]"]

    θ_M = Ω / (R * T.u) * (2 * Ω * E) / (9 * (1 - ν))
    c₀_cr = 0.0

    D_f = [params["$(Domain) diffusivity [m2.s-1]"](c[i], T.u) for i in 1:g.Nᵣ] # Diffusivities at cell centers
    Dₗ = [nothing; [D_face(D[i-1],D[i],Δr,Δr) for i in 2:g.Nᵣ]] # Left diffusivities
    Dᵣ = [[D_face(D[i], D[i+1],Δr,Δr) for i in 1:g.Nᵣ-1]; nothing] # Right diffusivities

    eqns = [
        # Diffusion with stress
        [σ[i] ~ 1 + θ_M * (c[i] - c₀_cr) for i in 1:g.Nᵣ]...
        [D[i] ~ params["$(Domain) diffusivity [m2.s-1]"](c[i], T.u)*σ[i] for i in 1:g.Nᵣ]...
        D_r ~ sum([D[i]*Vᵢ[i] for i in 1:g.Nᵣ])/sum(Vᵢ)
        
        c_avr ~ sum(c)/g.Nᵣ
        c_r ~ sum([c[i]*Vᵢ[i] for i in 1:g.Nᵣ])/sum(Vᵢ)
        c_surf ~ 1.5*c[end] - 0.5*c[end-1] # Surface concentration
        z ~ c_surf/params["Maximum concentration in $(domain) [mol.m-3]"] # Stoichiometry
        U₀ ~ params["$(Domain) OCP [V]"](z) # Open-circuit potential
        aₖ ~ 3*ϵₛ/params["$(Dmn) particle radius [m]"] # Specific surface area

        # Boundary condition center
        Dt(c[1]) ~ (Dᵣ[1]*Aᵣ[1]*(c[2] - c[1])/Δr)/Vᵢ[1]

        # Internal nodes
        [Dt(c[i]) ~ (Dᵣ[i]*Aᵣ[i]*(c[i+1] - c[i])/Δr - Dₗ[i]*Aₗ[i]*(c[i] - c[i-1])/Δr)/Vᵢ[i] for i in 2:g.Nᵣ-1]...

        # Boundary condition edge
        Dt(c[end]) ~ (-Aᵣ[end]*J.u/F - Dₗ[end]*Aₗ[end]*(c[end] - c[end-1])/Δr)/Vᵢ[end]

    ]


    # Event working
    events = [
        [
            c_surf ~ 0.01*params["Maximum concentration in $(domain) [mol.m-3]"],
            c_surf ~ 0.99*params["Maximum concentration in $(domain) [mol.m-3]"],
        ]=>(abort!,(;))
    ]

    System(eqns,t; name=name,systems=[J, T], continuous_events=events)
end