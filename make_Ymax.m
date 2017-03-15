function Ymax = make_Ymax(params, x0vec, inertial)

xk = x0vec(1); yk = x0vec(2); xdotk = x0vec(4); ydotk = x0vec(5);
phik = x0vec(8); thetak = x0vec(3);
rp = params.rp; rs = params.rs; rtol = params.rtol; gamma = params.gamma;
eta = params.eta; betaHIGH = params.betaHIGH; betaLOW = params.betaLOW;

if inertial==1
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
    
elseif inertial == 0
    Ymax = zeros(4,1);
    D = abs(xk-(rp+rs))+abs(yk);
    Ymax(1) = 0; Ymax(2) = 0;
    Ymax(3) = eta*D+betaHIGH;
    Ymax(4) = eta*D-betaLOW;
    
elseif inertial == 2
    Ymax = [0; 0];
    %Ymax = [0; 0; params.beta1 - eta*(rp+rs); -params.beta2 - eta*(rp+rs)];
    
elseif inertial == 3
    d = x0vec(2)/x0vec(1);
    Ymax = zeros(7,1);
    Ymax(1) = 0; Ymax(2) = 0;
    Ymax(3) = params.vnmax;
    Ymax(4) = params.vtmax;
    Ymax(5) = params.vtmax; 
    Ymax(6) = params.angtol+atan(d);
    Ymax(7) = params.angtol-atan(d);
end

