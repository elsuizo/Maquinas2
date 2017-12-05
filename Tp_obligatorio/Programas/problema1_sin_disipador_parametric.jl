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

step_power(t) = 3.871
pulse(t) = (175 * sin(2*π*100*t) + 175) / 2

family = @ode_def heat_simulation begin
   dT_c = (-T_c / C_c)*( 1 / (R_jc) + 1 / (R_ca))  + (T_j / (R_jc * C_c)) + (T_a / (R_ca * C_c))
   dT_j = (T_c/ (C_j * R_jc)) - T_j / (C_j * R_jc) + (step_power(t) / C_j)
end R_jc=>1.0 C_c=>6.8 R_ca=>30.0 C_j=>0.5 T_a=>30.0

T₀ = [0.0 0.0]
tspan = (0.0, 1500.0)
problems = [ODEProblem(heat_simulation(T_a=parameter), T₀, tspan) for parameter in 30:40]
solutions = solve.(problems, Vern7(), abstol=1/10^14, reltol=1/10^14)
for (num, solution) in enumerate(solutions)
   plot(solution.t, solution[2, :], label=latexstring("T_a="*"$(num + 29)"))
   legend(bbox_to_anchor=(1.09,1), loc="upper right", ncol=1)
   xlabel(L"t")
   ylabel(L"T")
end
