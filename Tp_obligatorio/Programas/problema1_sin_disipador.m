% definimos las constantes
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
q = 175;
f = @(t, T)[a_11 * T(1) + a_12 * T(2) + b; a_21 * T(1) + a_22 * T(2) - (q * C_j)];
[t, T_sol] = ode45(f, [0, 3600], [0, 0]);
plot(t,T_sol(:,1))

title('T(t)')
xlabel('t'), ylabel('y')
