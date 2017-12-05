#= -------------------------------------------------------------------------
# @file problema1_sin_disipador.jl
#
# @date 11/23/17 20:37:44
# @author Martin Noblia
# @email mnoblia@disroot.org
#
# @brief
# Simulacion de las ecuaciones difenciales del problema 1 sin disipador
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
# imports
using DifferentialEquations
#= using PyPlot =#
using Plots; pgfplots()
#=------------------------------------------------------------------------------
                           constants
------------------------------------------------------------------------------=#
const  C_j   =  0.5
const  C_c   =  6.8
const  R_jc  =  1.0
const  R_ca  =  30.0
const  T_a   =  30.0
const  T_j_max = 150.0
const  P_d  = (T_j_max - T_a) / (R_jc + R_ca)
a_11 = (-1 / C_c) * ((1 / R_jc) + (1 / R_ca))
a_12 = (1 / (R_jc * C_c))
a_21 = (1 / (C_j * R_jc))
a_22 = (-1/ (C_j * R_jc))
#=------------------------------------------------------------------------------
                        system matrix
------------------------------------------------------------------------------=#
A = [a_11 a_12
     a_21 a_22]

T₀ = [0.0, 0.0] # initial condition
tspan = (0.0, 1500.0) # simulation time
#=------------------------------------------------------------------------------
                           inputs
------------------------------------------------------------------------------=#
step_Pd(t) = (t >=0.0) ? P_d: 0.0
pulse_Pd(t) = (P_d * sin(2.0*π*100.0*t) + P_d) / 2.0
# function to modelate the system
f(t, T) = A * T + [0, pulse_Pd(t) / C_j] + [T_a / (R_ca * C_c), 0]
prob = ODEProblem(f, T₀, tspan)
sol = solve(prob)
final_value_T_c = sol[1, end]
final_value_T_j = sol[2, end]
println("el valor final es:", sol[2, end])
#=------------------------------------------------------------------------------
                           Ploting
------------------------------------------------------------------------------=#
plot(sol)
#= plot(sol.t, sol[1, :], label=L"T_c", linestyle="-.") =#
#= plot(sol.t, sol[2, :], label=L"T_j") =#
savefig("sample.tex")
#= axhline(final_value_T_j, color="k", linestyle="--") =#
#= text(0, final_value_T_j + 3, latexstring("T_J="*"$final_value_T_j"), fontsize=10) =#
#= xlabel(L"t\,[seg]") =#
#= ylabel(L"T\,[^{\circ}C]") =#
#= grid("on") =#
#= legend() =#
