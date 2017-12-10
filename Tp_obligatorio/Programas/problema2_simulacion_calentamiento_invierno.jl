#= -------------------------------------------------------------------------
# @file problema2_simulacion_calentamiento_invierno.jl
#
# @date 12/07/17 22:21:52
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
using DifferentialEquations
using PyPlot
#=------------------------------------------------------------------------------
                              constants
------------------------------------------------------------------------------=#
const  A      =  2.0
const  R      =  2.390057361376673e-5
const  cₑ     =  4.184e6
const  V      =  0.5
const  θ_ref  =  333.0
const  θ₀     =  281.0
const  θ_a    =  281.0
q = 1841012.0
K_p = 0.0001
G_inv = 0.4 * 836800.0
tspan = (0.0, 1000.0)
θ_a_f(t) = (100 < t > 300 ) ? θ_a + 100: θ_a
# calculo del q_comp
q_nec = ((θ_ref - θ_a) / (R)) - A * G_inv
# sin compensar en invierno
θ_dot(t, θ) = ((G_inv * A) / (cₑ * V)) + ((θ_a / (R * cₑ * V))) - ((θ / (R * cₑ * V)))
# compensado en invierno
θ_dot_comp(t, θ) = ((G_inv * A) / (cₑ * V)) + ((θ_a_f(t) / (R * cₑ * V))) - ((θ / (R * cₑ * V))) + (q_nec / (cₑ * V)) + (sol(t)-θ_ref) * K_p
θ_dot(t, θ) = ((q / (cₑ * V)) + (G_inv * A) / (cₑ * V)) + ((θ_a / (R * cₑ * V))) - ((θ / (R * cₑ * V)))
θ_dot(t, θ) = ((q / (cₑ * V)) + (G_inv * A) / (cₑ * V)) + ((θ_a / (R * cₑ * V))) - ((θ / (R * cₑ * V)))
prob = ODEProblem(θ_dot, θ₀, tspan)
prob_comp = ODEProblem(θ_dot_comp, θ₀, tspan)
sol = solve(prob, Tsit5(), reltol=1e-12,abstol=1e-12)
sol_comp = solve(prob_comp, Tsit5(), reltol=1e-12,abstol=1e-12)
#=------------------------------------------------------------------------------
                           plotting
------------------------------------------------------------------------------=#
plot(sol.t, sol[1, :] - 273.0, label=L"\theta_{inv}(t)", linewidth=3)
plot(sol_comp.t, sol_comp[1, :] - 273.0, label=L"\theta_{inv-comp}(t)", linewidth=3)
xlabel(L"t\,[seg]")
ylabel(L"T\,[^{\circ}C]")
grid("on")
matplotlib["rcParams"][:update](["font.size" => 14, "font.family" => "serif"])
legend()


