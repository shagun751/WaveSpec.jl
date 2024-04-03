#module run_wave


using   Revise
using   DrWatson
@quickactivate "WaveSpec"

using   WaveSpec
using   Plots
using   .Constants
using   .Jonswap
using   .WaveTimeSeries

pltFolder = datadir("output","res_")

h0 = 50 #still-water depth
H₀ = 7.11 #m
Tₚ = 12 #s
nω = 257

ω, S, A = jonswap(H₀, Tₚ; 
    plotflag=true, plotloc=pltFolder, nω = nω)

k = dispersionRelAng.(h0, ω; msg=false)
α = randomPhase(ω; seed=100)

sp = SpecStruct( h0, ω, S, A, k, α; Hs = H₀, Tp = Tₚ )

scatter(ω, k)
savefig(pltFolder*"dispersion.png")

scatter(ω, α)
savefig(pltFolder*"phase.png")

t=0.6
η, ϕ, u, w = waveAiry1D(h0, ω, A, k, α, t, 0.1, -0.1)
@show η, ϕ, u, w

t=0:0.1:1200
η, ϕ, u, w = waveAiry1D(h0, ω, A, k, α, t, 0.1, -0.1)

t=0.6
η, ϕ, u, w = waveAiry1D(sp, t, 0.1, -0.1)
@show η, ϕ, u, w
η, px, py = waveAiry1D_pPos(sp, t, 0.1, -0.1)
@show η, px, py

t=0:0.1:1200
η, ϕ, u, w = waveAiry1D(sp, t, 0.1, -0.1)

η, px, py = waveAiry1D_pPos(sp, t, 0.1, -0.1)

plot(t, η, dpi=330)
savefig(pltFolder*"ts_elevation.png")

plot(t, ϕ, dpi=330)
savefig(pltFolder*"ts_phi.png")

plot(t, u, dpi=330)
savefig(pltFolder*"ts_velu.png")

plot(t, w, dpi=330)
savefig(pltFolder*"ts_velw.png")

plot(t, px, dpi=330)
savefig(pltFolder*"ts_px.png")

plot(t, py, dpi=330)
savefig(pltFolder*"ts_py.png")


# end