#= -------------------------------------------------------------------------
# @file problema1_sin_disipador_parametric.jl
#
# @date 11/25/17 00:36:32
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

q(t) = 175

family = @ode_def heat_simulation begin
   dT_c = (-1 / C_c)*( 1 / R_jc + 1 / R_ca) * T_c + (1 / (R_jc * C_c)) * T_j + (T_a / R_ca * C_c)
   dT_j = (1 / (C_j * R_jc)) * T_c + (-1/ (C_j * R_jc)) * T_j - (q(t) * C_j)
end R_jc=>1.0 C_c=>6.8 R_ca=>30.0 C_j=>0.5 T_a=>30.0

T₀ = [0.0 0.0]
tspan = (0.0, 3600.0)
problems = [ODEProblem(heat_simulation(T_a=parameter), T₀, tspan) for parameter in 30:45]
solutions = solve.(problems)
for (num, solution) in enumerate(solutions)
   plot(solution[2, :], label=latexstring("T_a="*"$(num + 30)"))
   legend(bbox_to_anchor=(1.09,1), loc="upper right", ncol=1)
   xlabel(L"t")
   ylabel(L"T")
end
