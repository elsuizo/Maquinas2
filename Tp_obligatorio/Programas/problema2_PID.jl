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
import Sundials
using PyPlot
#=------------------------------------------------------------------------------
                              constants
------------------------------------------------------------------------------=#
const  A      =  2.0
const  R      =  2.390057361376673e-5
const  cₑ     =  4.184e6
const  V      =  0.5
const  θ_ref  =  333.0
const  θ₀     =  288.0
const  θₐ    =  293.0
G_inv = 0.4 * 836800.0
i(t) = 7
q(t) = (i(t))^2
function heat_panel(t, θ, dθ)
   dθ = q(t)/(cₑ * V) + (A * G_inv)/(cₑ * V * R) - (θ - θₐ) / (cₑ * V * R)
end

