function dTdt = system1(states)
   input = states(1);
   T_c   = states(2);
   T_j   = states(3);

   C_j   =  0.5;
   C_c   =  6.8;
   R_jc  =  1;
   R_ca  =  30;
   T_a   =  30;
   a_11 = (-1 / C_c) * ( 1 / R_jc + 1 / R_ca);
   a_12 = (1 / (R_jc * C_c));
   a_21 = (1 / (C_j * R_jc));
   a_22 = (-1/ (C_j * R_jc));
   b = (T_a / R_ca * C_c);
   dT_cdt = a_11 * T_c + a_12 * T_j + b;
   dT_jdt = a_21 * T_c + a_22 * T_j - (input * C_j);

   dTdt(1) = dT_cdt;
   dTdt(2) = dT_jdt;
end
