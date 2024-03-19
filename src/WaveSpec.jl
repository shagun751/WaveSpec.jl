module WaveSpec

using Revise 

include("Constants.jl")
include("Jonswap.jl")
include("Currents.jl")
include("WaveTimeSeries.jl")

export Constants, Jonswap, WaveTimeSeries
export Currents

end
