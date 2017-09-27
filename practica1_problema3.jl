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
using Plots; pyplot()

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
I₁(t) = (t ≥ 0) ? 5.2 : 0.0
I₂(t) = t
I₃(t) = sin(t)
#differential equation
f(t, T) = (I₁(t)^2 * Rₑ - π * D * h * (T - T_inf) - π * D * ϵ * σ *(T^4 - T_alr^4)) / (ρ * c *(π * (D^2)/4))
T₀ = 25.0 + 273.0
tsimu = (0.0, 240.0)
prob  = ODEProblem(f, T₀, tsimu)
solution = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

plot(solution.t, solution - 273)
