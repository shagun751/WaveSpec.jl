module WaveTimeSeries
# -------------------#
# This next section compiles the wave from 
# all seperate waves into one irregular wave

include("../src/constants.jl")

using .Constants

export waveAiry1D, waveStokes2nd1D

#---------- waveStokes2nd1D ----------#
function waveStokes2nd1D(h0::Real, ω::Real, A::Real, 
  k::Real, α::Real, t;
  x::Real=0, z::Real=0, x0::Real=0)

  αₘ = k*(x-x0) .- ω*t .+ α

  η = A * cos.(αₘ)

  ϕ = ifelse( ω > 0,
    g* A / ω * cosh(k*(h0+z)) / cosh(k*h0) * sin.(αₘ),
    0.0 * t )
  
  u = ifelse( ω > 0, 
    A * ω * cosh.(k*(h0+z)) / sinh(k*h0) * cos.(αₘ),
    0.0 * t )
  
  w = ifelse( ω > 0, 
    A * ω * sinh(k*(h0+z)) / sinh(k*h0) * sin.(αₘ),
    0.0 * t )

  return η, ϕ, u, w
  
end
#-------- End waveStokes2nd1D --------#


#------------- waveAiry1D ------------#
function waveAiry1D(h0::Real, ω, A, k, α, t;
  x::Real=0, z::Real=0, x0::Real=0)    

  params = (h0, ω, A, k, α, x0, x, z)
  res = waveAiry1D.( Ref(params), t)
  # Ref() is to ensure that broadcast operation is only
  # applied to `t`    

  if(res isa Array)
    η = [r for (r, nothing, nothing, nothing) in res]
    ϕ = [r for (nothing, r, nothing, nothing) in res]
    u = [r for (nothing, nothing, r, nothing) in res]
    w = [r for (nothing, nothing, nothing, r) in res]
  
  else
    η, ϕ, u, w = res
  end

  return η, ϕ, u, w
end


function waveAiry1D( params, t::Real)    

  h0, ω, A, k, α, x0, x, z = params

  αₘ = k*(x-x0) - ω*t + α

  η = sum( A .* cos.(αₘ) )

  ϕ = sum( ifelse.( ω .> 0,
    g* A ./ ω .* sin.(αₘ) .* cosh.(k*(h0+z)) ./ cosh.(k*h0),
    0.0 ) )
  
  u = sum( ifelse.( ω .> 0, 
    A .* ω .* cos.(αₘ) .* cosh.(k*(h0+z)) ./ sinh.(k*h0),
    0.0 ) )
  
  w = sum( ifelse.( ω .> 0, 
    A .* ω .* sin.(αₘ) .* sinh.(k*(h0+z)) ./ sinh.(k*h0),
    0.0 ) )

  return η, ϕ, u, w
end
#----------- End waveAiry1D ----------#


end 