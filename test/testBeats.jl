# Simulate single and two freq to observe the beats effect
using   WaveSpec
using   Plots
using   .Constants
using   .WaveTimeSeries

ω = [2*pi/9.9, 2*pi/10.1]
A = [1,1]
α = [0.0, 0.0]
k = dispersionRelAng.(30, ω)

t=0:0.1:500
η, ϕ, u, w = waveAiry1D(30, ω, A, k, α, t)
# η, ϕ, u, w = waveStokes2nd1D(30, ω[1], A[1], k[1], α[1], t)

@show k
