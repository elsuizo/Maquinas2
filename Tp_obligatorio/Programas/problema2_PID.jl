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
const  θ_ref  =  60.0 + 273.0
const θ₀_inicial     =  rand(273.0:290.0)
const  G_inv  = 0.4 * 836800.0
const  θₐ    =  293.0
q_nec = ((θ_ref - θₐ) / (R)) - A * G_inv
const  τ      = R * cₑ * V
α = cₑ * V
f(u) = u[2]
# funcion para generar una perturbacion aleatoria escalon
θ_a_f(t) = (100.0 < t < 300.0) ? θₐ - rand(1.0:10.0): θₐ

"""
Funcion que modela el problema 2
"""
function heat_panel(t, θ, dθ)
   dθ[1] = (f(u) / (α)) + ( ((A * G_inv)) / α - ((θ[1] - θₐ) / τ ) + q_nec / α)
end

function heat_panel_pert(t, θ, dθ)
   dθ[1] = (f(u) / (α)) + ( ((A * G_inv)) / α - ((θ[1] - θ_a_f(t)) / τ ) + q_nec / α)
end
#=------------------------------------------------------------------------------
                        constants of the algorithm
------------------------------------------------------------------------------=#
K_p = 1 / R
K_i = 837.0
tauI = 110
tauD = 0.0
t_sim = 500.0
δt  = .01
time_sim = linspace(0, 500, 10000)
set_point = zeros(time_sim) * θ_ref
error = zeros(time_sim) # zeros like t
θ = zeros(time_sim) # zeros like t
#= θ₀ = zeros(time_sim) # zeros like t =#
dθ = zeros(time_sim) # zeros like t
error_int = zeros(time_sim) # zeros like t
P = zeros(time_sim)
I = zeros(time_sim)
D = zeros(time_sim)
pv = zeros(time_sim)
controller_out = zeros(time_sim)
# Upper and Lower limits on OP
controller_max = 3e6
controller_min = 0.0

θ₀ = zeros(1)
θ₀[1] = θ₀_inicial
controller_min = 0
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
   controller_out[sample+1] = controller_out[1] + P[sample+1] + I[sample+1]
   if controller_out[sample+1] > controller_max  # check upper limit
      controller_out[sample+1] = controller_max
      error_int[sample+1] = error_int[sample+1] - error[sample+1] * Δt # anti-reset windup
   end
   if controller_out[sample+1] < controller_min  # check lower limit
      controller_out[sample+1] = controller_min
      error_int[sample+1] = error_int[sample+1] - error[sample+1] * Δt # anti-reset windup
   end
   ts = [time_sim[sample+1],time_sim[sample+2]]
   #= @show ts =#
   u[sample+2]= controller_out[sample+1]
   sol = Sundials.cvode(heat_panel_pert, θ₀,ts)
   #= sol = integrate.odeint(heat_panel, θ₀[1], ts) =#
   θ[sample+2] = sol[end]
   θ₀[1] = θ[sample+2]
   pv[sample+2] = θ[sample+2]
end

plot(time_sim[2:end], θ[2:end] - 273.0)
xlabel(L"t\,[seg]")
ylabel(L"T\,[^{\circ}C]")
grid("on")
matplotlib["rcParams"][:update](["font.size" => 14, "font.family" => "serif"])
text(0, θ₀_inicial - 273.0, latexstring("\\theta_{0}="*"$(θ₀_inicial-273.0)"), fontsize=10)
legend()
