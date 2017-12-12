#= -------------------------------------------------------------------------
# @file problema2_PID_simple.jl
#
# @date 12/09/17 17:10:03
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
using PyPlot

struct Panel_constants{T<:Number} <:Real
    A::T
    R::T
    cₑ::T
    V::T
    τ::T
    function Panel_constants{T}(A::T, R::T, cₑ::T, V::T) where T<:Number
       τ = R * cₑ * V
       new(A, R, cₑ, V, τ)
    end
end

Panel_constants(a::T, b::T, c::T, d::T) where T <: Number = Panel_constants{T}(a, b, c, d)

struct Problem_variables{T<:Number} <:Real
   θ₀::T
   θₐ::Function
   G::T
end

θ_ref = 60.0 + 273.0
θ_a = 293.0
G_inv = 0.4 * 836800.0
θ_a_f(t) = 293.0
constants = Panel_constants(2.0, 2.390e-5, 4.184e6, 0.5)
problem_vars = Problem_variables(288.0,θ_a_f, G_inv)
q_nec = ((θ_ref - θ_a) / (constants.R)) - constants.A * G_inv
#------------------------------------------------------------------------------
#                            function
#------------------------------------------------------------------------------
function heat_panel(t, constants::Panel_constants, problem::Problem_variables)
   problem_vars.θ₀ * exp(-t / constants.τ) + ((problem_vars.θₐ(t) + constants.R*constants.A*problem_vars.G) * (1 - exp(-t / constants.τ)))
end
K_p = 0.001
t = 0:0.01:500
sample_past = 0.0
θ = [heat_panel(sample, constants, problem_vars) + q_nec  for sample in t]
plot(t, θ-273)

