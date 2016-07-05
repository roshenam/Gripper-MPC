function Ymax = make_Ymax(params, x0vec)

xk = x0vec(1); yk = x0vec(2); xdotk = x0vec(4); ydotk = x0vec(5); 
phik = x0vec(8); thetak = x0vec(3); 
rp = params.rp; rs = params.rs; rtol = params.rtol; gamma = params.gamma;
eta = params.eta; betaHIGH = params.betaHIGH; betaLOW = params.betaLOW;
Ymax = zeros(4,1);
%Ymax = zeros(2,1);
L1 = sin(phik+gamma)/((rp-rtol)*sin(gamma));
L2 = cos(phik+gamma)/((rp-rtol)*sin(gamma));
L3 = sin(phik-gamma)/((rp-rtol)*sin(gamma));
L4 = cos(phik-gamma)/((rp-rtol)*sin(gamma));
D = abs(xk - (rp+rs)*cos(phik)) + abs(yk - (rp+rs)*sin(phik));
%update = (xdotk*sin(thetak) + ydotk*cos(thetak))*thetak;
Ymax(1) = 1 + phik*(L2*xk + L1*yk);
Ymax(2) = 1 + phik*(-L4*xk - L3*yk);
Ymax(3) = eta*D + betaHIGH;
Ymax(4) = eta*D - betaLOW;

