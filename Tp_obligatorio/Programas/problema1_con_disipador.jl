#= -------------------------------------------------------------------------
# @file problema1_con_disipador.jl
#
# @date 12/05/17 19:43:54
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
# problema con disipador
#=------------------------------------------------------------------------------
                        constants
------------------------------------------------------------------------------=#
const  C_j   =  0.5
const  C_c   =  6.8
const  C_s   =  39.6
const  R_jc  =  1.0
const  R_ca  =  30.0
const  R_sa  =  4.5
const  R_cs  =  0.4
const  T_a   =  30.0
const  T_j_max = 150.0
#= const  P_d  = (T_j_max - T_a) / (R_jc + (R_ca / (R_cs + R_sa))) =#
P_d = 23.02

a_11 = (-1 / C_c) * ((1 / R_jc) + (1 / R_ca) + (1 / R_cs))
a_12 = (1 / (R_jc * C_c))
a_13 = (1 / (R_cs * C_c))

a_21 = (1 / (C_j * R_jc))
a_22 = (-1/ (C_j * R_jc))
a_23 = 0.0

a_31 = (1 / (R_cs * C_s))
a_32 = 0.0
a_33 = (-1 / C_s) * ((1 / R_cs) + (1 / R_sa))

A = [a_11 a_12 a_13
     a_21 a_22 a_23
     a_31 a_32 a_33]

T₀ = [0.0, 0.0, 0.0] # initial condition
tspan = (0.0, 1500.0) # simulation time
#=------------------------------------------------------------------------------
                           inputs
------------------------------------------------------------------------------=#
step_Pd(t) = (t >=0.0) ? P_d: 0.0
pulse_Pd(t) = (P_d * sin(2.0*π*100.0*t) + P_d) / 2.0
# function to modelate the system
f(t, T) = A * T + [0, step_Pd(t) / C_s, 0] + [T_a / (R_ca * C_c), 0, T_a / (R_sa * C_s)]
prob = ODEProblem(f, T₀, tspan)
sol = solve(prob)
final_value_T_c = sol[1, end]
final_value_T_j = sol[2, end]
println("el valor final de T_J es:", sol[2, end])
#=------------------------------------------------------------------------------
                           Ploting
------------------------------------------------------------------------------=#
plot(sol.t, sol[1, :], label=L"T_c", linewidth=3)
plot(sol.t, sol[2, :], label=L"T_j", linewidth=2)
plot(sol.t, sol[3, :], label=L"T_s", linewidth=2)
axhline(final_value_T_j, color="k", linestyle="--")
text(0, final_value_T_j + 3, latexstring("T_J="*"$final_value_T_j"), fontsize=10)
xlabel(L"t\,[seg]")
ylabel(L"T\,[^{\circ}C]")
grid("on")
matplotlib["rcParams"][:update](["font.size" => 14, "font.family" => "serif"])
legend()
