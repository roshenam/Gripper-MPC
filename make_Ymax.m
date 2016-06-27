function Ymax = make_Ymax(params, xk, yk, phik)

rp = params.rp; rtol = params.rtol; gamma = params.gamma;
Ymax = zeros(2,1);
L1 = sin(phik+gamma)/((rp-rtol)*sin(gamma));
L2 = cos(phik+gamma)/((rp-rtol)*sin(gamma));
L3 = sin(phik-gamma)/((rp-rtol)*sin(gamma));
L4 = cos(phik-gamma)/((rp-rtol)*sin(gamma));
Ymax(1) = 1 + phik*(L2*xk + L1*yk);
Ymax(2) = 1 + phik*(-L4*xk - L3*yk);
