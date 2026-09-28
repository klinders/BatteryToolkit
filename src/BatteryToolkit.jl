module BatteryToolkit

using ModelingToolkit
using ModelingToolkitStandardLibrary
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkitStandardLibrary.Electrical

include("ParameterSets/Base.jl")
include("helpers.jl")
include("ParameterSets/Chen2020.jl")
include("ParameterSets/OKane2022.jl")
include("ParameterSets/Ruihe2025.jl")

export Chen2020, OKane2022, Ruihe2025, BatteryParameters

include("CellModels/SPMe/SPMe.jl")

export SPMe

include("ThermalModels/ConstantTemperature.jl")
include("ThermalModels/AmbientTemperature.jl")

export ConstantTemperature

include("experiment.jl")

export Experiment,AbstractStep,PowerStep,CurrentStep,DriveStep,RestStep,ChargeStep,step!, get_p0

include("solvers.jl")
include("PackModels/SingleCellPack.jl")
include("PackModels/MultiCellPack.jl")

export SingleCellPack, MultiCellPack, simulate

end  # module BatteryToolkit