#= -------------------------------------------------------------------------
# @file problema2_PID.jl
#
# @date 12/08/17 01:24:52
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
using Sundials
using PyPlot
#=------------------------------------------------------------------------------
                              constants
------------------------------------------------------------------------------=#
const  A      =  2.0
const  R      =  2.390057361376673e-5
const  cₑ     =  4.184e6
const  V      =  0.5
const  θ_ref  =  333.0
const  θ₀     =  288.0
const  θₐ    =  293.0
G_inv = 0.4 * 836800.0


function heat_panel(t, θ, dθ)
   dθ[1] = u[1] /(cₑ * V) + (A * G_inv)/(cₑ * V * R) - (θ[1] - θₐ) / (cₑ * V * R)
end
#=------------------------------------------------------------------------------
                        constants of the algorithm
------------------------------------------------------------------------------=#
K_c = 1
tauI = 0.913444964569
tauD = 0.0
t_sim = 100
δt  = .1
time_sim = linspace(0, t_sim)
set_point = zeros(time_sim) # zeros like t
error = zeros(time_sim) # zeros like t
θ = zeros(time_sim) # zeros like t
dθ = zeros(time_sim) # zeros like t
error_int = zeros(time_sim) # zeros like t
P = zeros(time_sim)
I = zeros(time_sim)
D = zeros(time_sim)
pv = zeros(time_sim)
controller_out = zeros(time_sim)
# Upper and Lower limits on OP
controller_max = 350.0
controller_min = 250.0

u = ones(time_sim) * θ_ref
for sample in 0:t_sim
   Δt = time_sim[sample+1] - time_sim[sample+2]
   error[sample+1] = set_point[sample+1] - pv[sample+1]
   P[sample+1] = K_c * error[sample+1]
   I[sample+1] = K_c/tauI * error_int[sample+1]
   D[sample+1] = - K_c * tauD * dθ[sample+1]
   controller_out[sample+1] = controller_out[1] + P[sample+1] + I[sample+1] + D[sample+1]
   if controller_out[sample+1] > controller_max  # check upper limit
      controller_out[sample+1] = controller_max
      error_int[sample+1] = error_int[sample+1] - error[sample+1] * Δt # anti-reset windup
   end
   if controller_out[sample+1] < controller_min  # check lower limit
      controller_out[sample+1] = controller_min
      error_int[sample+1] = error_int[sample+1] - error[sample+1] * Δt # anti-reset windup
   end
   ts = [time_sim[sample+1],time_sim[sample+2]]
   u[sample+2]= controller_out[sample+1]
   sol = Sundials.cvode(heat_panel, θ₀,ts)
   θ[sample+2] = sol[end]
   θ₀ = θ[sample+2]
   pv[sample+2] = θ[sample+2]
end
#= controller_out[len(t)-1] = controller_out[len(t)-2] =#
#= ie[len(t)-1] = ie[len(t)-2] =#
#= P[len(t)-1] = P[len(t)-2] =#
#= I[len(t)-1] = I[len(t)-2] =#
#= D[len(t)-1] = D[len(t)-2] =#
#=  =#

