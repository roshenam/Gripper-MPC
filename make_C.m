function C = make_C(params, x0vec, slack)

rp = params.rp; rtol = params.rtol; gamma = params.gamma;
xk = x0vec(1); yk = x0vec(2); xdotk = x0vec(3); ydotk = x0vec(4); 
phik = x0vec(6); thetak = x0vec(7); 
if slack
    C = zeros(4,8);
    L1 = sin(phik+gamma)/((rp-rtol)*sin(gamma));
    L2 = cos(phik+gamma)/((rp-rtol)*sin(gamma));
    L3 = sin(phik-gamma)/((rp-rtol)*sin(gamma));
    L4 = cos(phik-gamma)/((rp-rtol)*sin(gamma));
    
    % Cone constraints
    C(1,1) = L1;
    C(1,2) = -L2;
    C(2,1) = -L3;
    C(2,2) = L4;
    C(1,6) = L2*xk+L1*yk;
    C(2,6) = -L4*xk-L3*yk;
    
    % Normal velocity constraints
    C(3,3) = cos(thetak); 
    C(3,4) = sin(thetak); 
    C(4,3) = -C(3,3);
    C(4,4) = -C(3,4);

    
else
    C = zeros(2,10);
    L1 = sin(phik+gamma)/((rp-rtol)*sin(gamma));
    L2 = cos(phik+gamma)/((rp-rtol)*sin(gamma));
    L3 = sin(phik-gamma)/((rp-rtol)*sin(gamma));
    L4 = cos(phik-gamma)/((rp-rtol)*sin(gamma));
    
    C(1,1) = L1;
    C(1,2) = -L2;
    C(2,1) = -L3;
    C(2,2) = L4;
end

    