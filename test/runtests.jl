using WaveSpec
using Test

@time @testset "JONSWAP" begin include("testJonswap.jl") end
@time @testset "Beats" begin include("testBeats.jl") end
