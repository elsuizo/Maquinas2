#= -------------------------------------------------------------------------
# @file simulation_PID.jl
#
# @date 12/12/17 19:19:33
# @author Martin Noblia
# @email mnoblia@disroot.org
#
# @brief
#
# @detail
#
#  Licence:
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License

---------------------------------------------------------------------------=#
module SimulationPID

#------------------------------------------------------------------------------
#                            parameters
#------------------------------------------------------------------------------
const  A      =  2.0
const  R      =  2.390057361376673e-5
const  cₑ     =  4.184e6
const  V      =  0.5
const  θ_ref  =  333.0
const  θ₀     =  281.0
const  θ_a    =  281.0
const  G_inv  = 0.4 * 836800.0
const  α      = cₑ * V
const  τ      = R * cₑ * V
const  q_nec = ((θ_ref - θ_a) / (R)) - A * G_inv

#------------------------------------------------------------------------------
#                            function system
#------------------------------------------------------------------------------


function update()
end

function update(system::Function, initial_condition::Number)
  @assert stop > start
  # Pre-allocate the timeseries matrix
  t = (start, stop)
  t_keep = collect(start:1.0:stop)

  # Pre-assign function
  θ_dot_comp(t, θ) = (G_inv * A / α) + (θ_a / τ) - (θ / τ) + (q_nec / α)

  # Perform the actual integration
  prob = ODEProblem(f, biomass, t)
  sol = solve(prob, saveat=t_keep, dense=false, save_timeseries=false)

end

function simulate(system::Function, start::Real, stop::Real, δ::Number)
   @assert stop > start
   output = allocate_simulation_out()
   while N > 1
      for time in start:δ:stop
         
      end
   end
end

function simulate(p, biomass; start::Int64=0, stop::Int64=500, use::Symbol=:stiff)
  @assert length(biomass) == size(p[:A],1)
  @assert use ∈ vec([:stiff :nonstiff])

  # Pre-allocate the timeseries matrix
  t = (float(start), float(stop))
  t_keep = collect(start:1.0:stop)

  # Pre-assign function
  f(t, y) = dBdt(t, y, p)

  # Perform the actual integration
  prob = ODEProblem(f, biomass, t)
  sol = solve(prob, saveat=t_keep, dense=false, save_timeseries=false, alg_hints=[use])

  output = Dict{Symbol,Any}(
  :p => p,
  :t => sol.t,
  :B => hcat(sol.u...)'
  )

  return output

end

end # end of SimulationPID
