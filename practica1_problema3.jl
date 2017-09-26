#= -------------------------------------------------------------------------
# @file practica1_problema3.jl
#
# @date 09/22/17 13:41:52
# @author Martin Noblia
# @email martin.noblia@openmailbox.org
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
using Plots

#---------------------------------------------------------------------------
# constants
#---------------------------------------------------------------------------
const ϵ = 0.8
const Rₑ = 0.4
const D = 1e-3
const T_alr = 300
const T_inf = 300
const h = 100
const σ = 5.67e-8
const ρ = 8933
const c = 385
const I = 3
#differential equation
f(t, T) = (I^2 * Rₑ - π * D * h * (T - T_inf) - π * D * ϵ * σ *(T^4 - T_alr^4)) / (ρ * c *(π * (D^2)/4))
T₀ = 0.0
tspan = (0.0, 3.0)
prob  = ODEProblem(f, T₀, tspan)
solution = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

#= plot(solution,linewidth=5,title="Solution to the linear ODE with a thick line", =#
#=      xaxis="Time (t)",yaxis="u(t) (in μm)",label="My Thick Line!") # legend=false =#
