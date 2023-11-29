# Generate JONSWAP spectrum

using   WaveSpec
using   Plots
using   .Constants
using   .Jonswap
using   .WaveTimeSeries

h0 = 50 #still-water depth
H₀ = 7.11 #m
Tₚ = 12 #s

ω, S, A = jonswap(H₀, Tₚ)

k = dispersionRelAng.(h0, ω; msg=false)
α = randomPhase(ω; seed=100)

t=0.6
η, ϕ, u, w = waveAiry1D(h0, ω, A, k, α, t)

t=0:0.1:0.6
η, ϕ, u, w = waveAiry1D(h0, ω, A, k, α, t)

