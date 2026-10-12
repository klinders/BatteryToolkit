
"""
    ElectrolyteParameters

Parameters describing lithium-ion battery electrolyte properties.

# Fields
- `Lₚ::Float64`: Positive electrode thickness (m)
- `Lₛ::Float64`: Separator thickness (m)
- `Lₙ::Float64`: Negative electrode thickness (m)
- `Dₑ::Function`: Electrolyte diffusivity as function of concentration (m²/s)
- `σₑ::Function`: Electrolyte conductivity as function of concentration (S/m)
- `c₀::Float64`: Initial electrolyte concentration (mol/m³)
- `cₜ::Float64`: Typical electrolyte concentration (mol/m³)
- `t₊::Float64`: Transference number (cation fraction in current)
- `ϵₚ::Float64`: Positive electrode porosity
- `ϵₛ::Float64`: Separator porosity
- `ϵₙ::Float64`: Negative electrode porosity
- `bₚ::Float64`: Positive electrode Bruggeman coefficient
- `bₛ::Float64`: Separator Bruggeman coefficient
- `bₙ::Float64`: Negative electrode Bruggeman coefficient
"""
Base.@kwdef mutable struct ElectrolyteParameters
    Lₚ # Length of the electrode in m
    Lₛ # Length of the separator in m
    Lₙ # Length of the negative electrode in m
    Dₑ # Diffusivity in m^2*s^-1
    σₑ # Conductivity in S*m^-1
    c₀ # Initial electrolyte concentration in mol*m^-3
    cₜ # Typical Electrolyte concentration in mol*m^-3
    t₊ # transer number
    ϵₚ # Porosity of the positive electrode
    ϵₛ # Porosity of the separator
    ϵₙ # Porosity of the negative electrode
    bₙ # Negative electrode Bruggeman coefficient
    bₛ # Separator Bruggeman coefficient
    bₚ # Positive electrode Bruggeman coefficient
end

"""
    SideReactionParameters

Parameters for secondary reactions (SEI growth, lithium plating) on electrode surfaces.

# Fields
- `name::Symbol`: Reaction identifier (e.g., `:sei`)
- `k::Float64`: Reaction rate coefficient (A/m²)
- `α::Float64`: Transfer coefficient (charge transfer kinetics parameter)
- `M::Float64`: Molar mass of reaction product (kg/mol)
- `z::Int`: Number of electrons transferred per product molecule
- `ρ::Float64`: Product film density (kg/m³)
- `σ::Float64`: Product film conductivity (S/m)
- `U::Float64`: Open circuit potential of reaction (V)
- `Lf₀::Float64`: Initial film thickness (m)
- `V̄::Float64`: Partial molar volume of product (m³/mol)
- `R::Float64`: Film resistivity (Ω·m)
- `j_sei₀::Float64`: Exchange current density (A/m²)
- `c::Function`: Concentration dependence function
- `D_sol::Float64`: Solvent diffusivity in film (m²/s)
- `c_sol::Float64`: Solvent concentration (mol/m³)
"""
Base.@kwdef mutable struct SideReactionParameters
    name            # Name of the side reaction
    k               # Reaction rate
    α               # Side reaction transfer coefficient
    M               # Molar mass of SR product
    z               # Ratio of Li to SEI moles
    ρ               # Density of SR product
    σ               # Conductivity in the SEI layer
    U               # Open circuit potential of SR
    Lf₀             # Initial thickness of SR film
    V̄               # SEI Partial molar volume
    R               # SEI Resistivity
    j_sei₀          # Reaction exchange current
    c::Function     # Concentration dependence function
    D_sol           # Solvent diffusivity in m^2/s
    c_sol           # Solvent concentration in mol/m^3
    E_sei           # Activation energy for SEI growth in J/mol
    T_ref           # Reference temperature in K
    α_plating       # Li plating transfer coefficient
    k_plating       # Li plating reaction rate
    V_plating       # Partial molar volume of lithium [m3.mol-1]
end

"""
    SolidParticleParameters

Parameters describing lithium-ion electrode (active material particle) properties.

# Fields
- `Rₖ::Float64`: Particle radius (m)
- `aₖ::Float64`: Interfacial area per unit volume (m⁻¹)
- `Dₖ::Function`: Solid-state diffusivity as function of concentration (m²/s)
- `σₖ::Float64`: Electronic conductivity (S/m)
- `c₀::Float64`: Initial lithium concentration (mol/m³)
- `c₊::Float64`: Maximum lithium concentration (mol/m³)
- `Uₖ::Function`: Open circuit potential as function of stoichiometry (V)
- `j0::Function`: Exchange current density as function of electrolyte concentration, surface concentration, maximum concentration, and temperature (A/m²)
- `L_sei₀::Float64`: Initial SEI thickness (m)
- `side_reactions::Vector{SideReactionParameters}`: Secondary reactions on this electrode (default: [])
"""
Base.@kwdef mutable struct SolidParticleParameters
    Rₖ # Radius of the electrode in m
    ϵₛ # Active material volume fraction
    aₖ # Surface area density in m^-1
    Dₖ # electrode diffusivity in m^2*s^-1
    σₖ # Conductivity in S*m^-1
    c₀ # Initial electrode concentration in mol*m^-3
    c₊ # Maximum electrode concentration in mol*m^-3
    Uₖ # Open-circuit potential in V
    j0 # Exchange current density function
    L_sei₀ # Initial thickness of SEI film
    ρ_cr # Crack density
    w_cr # Crack width
    Ω # parial molar volume
    E # Youngs modulus
    ν # poissons ratio
    β_LAM # LAM rate constant
    m_LAM # LAM stress exponent
    stress_critical # Critical stress for LAM
    
    side_reactions::Vector{SideReactionParameters} = SideReactionParameters[]
end

"""
    BatteryParameters

Complete parameter set for a lithium-ion battery cell.

Aggregates all electrochemical and geometric parameters needed for SPMe simulations,
including positive/negative electrode, electrolyte, and pack configuration.

# Fields
- `p::SolidParticleParameters`: Positive electrode (cathode) parameters
- `n::SolidParticleParameters`: Negative electrode (anode) parameters
- `e::ElectrolyteParameters`: Electrolyte parameters
- `Hcc::Float64`: Current collector height (m)
- `Wcc::Float64`: Current collector width (m)
- `n_el::Int`: Number of parallel electrode pairs
- `Q₀::Float64`: Nominal cell capacity (Ah)
- `Vmin::Float64`: Minimum safe voltage (V)
- `Vmax::Float64`: Maximum safe voltage (V)

# Example
```julia
params = Chen2020()  # Pre-configured parameter set
sys = SPMe(params=params)
```
"""
Base.@kwdef mutable struct BatteryParameters
    p::SolidParticleParameters # Parameters for the positive electrode
    n::SolidParticleParameters # Parameters for the negative electrode
    e::ElectrolyteParameters # Parameters for the electrolyte
    
    Hcc # Current collector height
    Wcc # Current collector width
    n_el # Number of parallel electrodes
    Q₀ # Original battery capacity in Ah

    Vmin
    Vmax

    sei
    plating
    cracking_n
    cracking_p
    lam_n
    lam_p
    dissolution
end