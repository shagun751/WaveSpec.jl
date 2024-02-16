module WaveTimeSeries
# -------------------#
# This next section compiles the wave from 
# all seperate waves into one irregular wave

using Revise
# using DrWatson
# @quickactivate "WaveSpec"

using WaveSpec.Constants

export waveAiry1D


#------------- waveAiry1D ------------#

function waveAiry1D(sp::SpecStruct, t, x, z; 
  x0::Real = 0)    

  return waveAiry1D(sp.h0, sp.ω, sp.A, sp.k, sp.α, 
    t, x, z; 
    x0 = x0)
  
end


function waveAiry1D(h0::Real, ω, A, k, α, 
  t::AbstractArray{<:AbstractFloat},
  x, z; x0::Real=0)    

  params = (h0, ω, A, k, α, x0, x, z)
  res = waveAiry1D.( Ref(params), t)
  # Ref() is to ensure that broadcast operation is only
  # applied to `t`    

  η = [r for (r, nothing, nothing, nothing) in res]
  ϕ = [r for (nothing, r, nothing, nothing) in res]
  u = [r for (nothing, nothing, r, nothing) in res]
  w = [r for (nothing, nothing, nothing, r) in res]    

  return η, ϕ, u, w
end


function waveAiry1D(h0::Real, ω, A, k, α, t::Real,
  x, z; x0::Real=0)    

  params = (h0, ω, A, k, α, x0, x, z)
  res = waveAiry1D( params, t)  
  
  η, ϕ, u, w = res

  return η, ϕ, u, w
end


# function waveAiry1D(h0::Real, ω, A, k, α, t,
#   x, z; x0::Real=0)    

#   params = (h0, ω, A, k, α, x0, x, z)
#   res = waveAiry1D.( Ref(params), t)
#   # Ref() is to ensure that broadcast operation is only
#   # applied to `t`    

#   if(res isa Array)
#     η = [r for (r, nothing, nothing, nothing) in res]
#     ϕ = [r for (nothing, r, nothing, nothing) in res]
#     u = [r for (nothing, nothing, r, nothing) in res]
#     w = [r for (nothing, nothing, nothing, r) in res]
  
#   else
#     η, ϕ, u, w = res
#   end

#   return η, ϕ, u, w
# end


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