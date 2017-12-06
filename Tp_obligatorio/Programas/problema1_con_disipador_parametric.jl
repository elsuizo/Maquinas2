#= -------------------------------------------------------------------------
# @file problema1_con_disipador_parametric.jl
#
# @date 12/06/17 16:18:05
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
# imports
using DifferentialEquations
using PyPlot
#=------------------------------------------------------------------------------
                  sistem parametric equations
------------------------------------------------------------------------------=#
family = @ode_def heat_simulation begin
   dT_c = (-T_c / C_c)*( 1 / (R_ca) + 1 / (R_jc) + 1 / (R_cs))  + (T_j / (R_jc * C_c)) + (T_s / (R_cs * C_c))+ (T_a / (R_ca * C_c))
   dT_j = (T_c/ (R_jc * C_j)) - (T_j / (C_j * R_jc)) + (step_power(t) / C_j)
   dT_s = (T_c / (R_cs * C_s)) - (T_s / C_s) * (1 / (R_cs) + 1 / (R_sa)) + (T_a / (R_sa * C_s))
end R_jc=>1.0 C_c=>6.8 R_ca=>30.0 C_j=>0.5 T_a=>30.0 C_s=>39.6 R_cs=>0.4 R_sa=>4.5


