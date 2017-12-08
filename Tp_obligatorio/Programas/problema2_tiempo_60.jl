#= -------------------------------------------------------------------------
# @file problema2_tiempo_60.jl
#
# @date 12/07/17 18:00:56
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
#=------------------------------------------------------------------------------
                              constants
------------------------------------------------------------------------------=#
const  A      =  2.0
const  R      =  2.390057361376673e-5
const  cₑ     =  4.184e6
const  V      =  0.5
const  θ_ref  =  333.0
const  θ₀     =  288.0
const  θ_a    =  293.0
G = (θ_ref - θ_a) / (R * A)
tspan = (0.0, 1000.0)
θ_dot(t, θ) = ((G * A) / (cₑ * V)) + ((θ_a / (R * cₑ * V))) - ((θ / (R * cₑ * V)))
prob = ODEProblem(θ_dot, θ₀, tspan)
condition(t,θ,integrator) = (θ - 333.0)
affect!(integrator) = terminate!(integrator)
cb = ContinuousCallback(condition,affect!)
sol = solve(prob, callback=cb)
