#= -------------------------------------------------------------------------
# @file problema1_sin_disipador_sundials.jl
#
# @date 11/27/17 18:07:28
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
import Sundials
using PyPlot

const  C_j   =  0.5
const  C_c   =  6.8
const  R_jc  =  1.0
const  R_ca  =  30
const  T_a   =  30
P_d  = (150 - 30) / (1+30)
a_11 = (-1 / C_c) * ( ( 1 / R_jc ) + ( 1 / R_ca ))
a_12 = (1 / (R_jc * C_c))
a_21 = (1 / (C_j * R_jc))
a_22 = (-1/ (C_j * R_jc))

step(t) = (t >=0.0) ? P_d: 0.0
pulse(t) = (P_d * sin(2.0*π*100.0*t) + P_d) / 2.0
"""
docs
"""
function system(t, T, dT)
   dT[1] = a_11 * T[1] + a_12 * T[2] + (T_a / (R_ca * C_c))
   dT[2] = a_21 * T[1] + a_22 * T[2] + (pulse(t) / C_j)
end

T₀ = [0.0, 0.0] # initial condition
t = linspace(0, 1300, 1000000) # time vector
sol = Sundials.cvode(system, T₀, collect(t))
# plot the results
plot(t, sol[:,1], label=L"T_c")
plot(t, sol[:,2], label=L"T_j")
xlabel(L"t")
ylabel(L"T")
legend()
