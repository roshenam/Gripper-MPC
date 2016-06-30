function C = make_C(params, xk, yk, phik, slack)

rp = params.rp; rtol = params.rtol; gamma = params.gamma;

if slack
    C = zeros(2,8);
    L1 = sin(phik+gamma)/((rp-rtol)*sin(gamma));
    L2 = cos(phik+gamma)/((rp-rtol)*sin(gamma));
    L3 = sin(phik-gamma)/((rp-rtol)*sin(gamma));
    L4 = cos(phik-gamma)/((rp-rtol)*sin(gamma));
    
    C(1,1) = L1;
    C(1,2) = -L2;
    C(2,1) = -L3;
    C(2,2) = L4;
    C(1,6) = L2*xk+L1*yk;
    C(2,6) = -L4*xk-L3*yk;
    
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

    