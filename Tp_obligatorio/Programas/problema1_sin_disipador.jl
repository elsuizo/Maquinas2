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
#= using Plots =#
using PyPlot
# defining the constants
const  C_j   =  0.5
const  C_c   =  6.8
const  R_jc  =  1
const  R_ca  =  30
const  T_a   =  30

a_11 = (-1 / C_c)*( 1 / R_jc + 1 / R_ca)
a_12 = (1 / (R_jc * C_c))
a_21 = (1 / (C_j * R_jc))
a_22 = (-1/ (C_j * R_jc))
A = [a_11 a_12
     a_21 a_22]

T₀ = [0.0, 0.0]
tspan = (0.0, 3600.0)
q(t) = 175
f(t, T) = A * T - [0, q(t) * C_j] + [(1 / R_ca * C_c) * T_a, 0]
prob = ODEProblem(f, T₀, tspan)
sol = solve(prob)
plot(sol.t, sol[1, :])
