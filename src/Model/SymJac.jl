syms rho_p rho_m
syms E_p E_m
syms K_p K_m
syms Phi_p Phi_m
syms dz
p_p=E_p-rho_p*(K_p+Phi_p);
p_m=E_m-rho_m*(K_m+Phi_m);
D=2/(rho_p+rho_m)*(p_p-p_m)/dz;

dDdrho_p=simplify(diff(D,rho_p));
dDdrho_m=simplify(diff(D,rho_m));

aa=3;