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
θ₀     =  [288.0]
const  θₐ    =  293.0
G_inv = 0.4 * 836800.0

f(u) = u[2]
function heat_panel(t, θ, dθ)
   dθ[1] = (f(u) / (cₑ * V)) + ((A * G_inv))/((cₑ * V * R)) - ((θ[1] - θₐ)) / ((cₑ * V * R))
end
#=------------------------------------------------------------------------------
                        constants of the algorithm
------------------------------------------------------------------------------=#
K_p = 41840
K_i = 836.8
tauI = 110
tauD = 0.0
t_sim = 1000
δt  = .01
time_sim = collect(0:δt:t_sim)
set_point = ones(time_sim) * θ_ref
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
controller_max = 600.0
controller_min = -100.0

pv[1] = θ_ref
u = ones(time_sim) * θ_ref
for sample in 0:length(time_sim) - 2
   Δt = time_sim[sample+1] - time_sim[sample+2]
   error[sample+1] = set_point[sample+1] - pv[sample+1]
   if sample >= 1  # calculate starting on second cycle
     error_int[sample+1] = error_int[sample] + error[sample+1] * Δt
   end
   P[sample+1] = K_p * error[sample+1]
   I[sample+1] = K_i* error_int[sample+1]
   controller_out[sample+1] = controller_out[1] + P[sample+1] + I[sample+1] + D[sample+1]
   #= if controller_out[sample+1] > controller_max  # check upper limit =#
   #=    controller_out[sample+1] = controller_max =#
   #=    error_int[sample+1] = error_int[sample+1] - error[sample+1] * Δt # anti-reset windup =#
   #= end =#
   #= if controller_out[sample+1] < controller_min  # check lower limit =#
   #=    controller_out[sample+1] = controller_min =#
   #=    error_int[sample+1] = error_int[sample+1] - error[sample+1] * Δt # anti-reset windup =#
   #= end =#
   ts = [time_sim[sample+1],time_sim[sample+2]]
   #= @show ts =#
   u[sample+2]= controller_out[sample+1]
   sol = Sundials.cvode(heat_panel, θ₀,ts)
   @show sol
   θ[sample+2] = sol[end]
   θ₀[1] = θ[sample+2]
   pv[sample+2] = θ[sample+2]
end

#= plot() =#
