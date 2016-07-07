function C = make_C(params, x0vec, slack)

rp = params.rp; rtol = params.rtol; gamma = params.gamma;
xk = x0vec(1); yk = x0vec(2);
phik = x0vec(8);
if slack
    C = zeros(4,8);
    %C = zeros(2,8);
    L1 = sin(phik+gamma)/((rp-rtol)*sin(gamma));
    L2 = cos(phik+gamma)/((rp-rtol)*sin(gamma));
    L3 = sin(phik-gamma)/((rp-rtol)*sin(gamma));
    L4 = cos(phik-gamma)/((rp-rtol)*sin(gamma));
    
    % Cone constraints
    C(1,1) = L1;
    C(1,2) = -L2;
    C(2,1) = -L3;
    C(2,2) = L4;
    C(1,8) = L2*xk+L1*yk;
    C(2,8) = -L4*xk-L3*yk;
    
    % Normal velocity constraints
    C(3,4) = -cos(phik);
    C(3,5) = -sin(phik);
    C(4,4) = -C(3,4);
    C(4,5) = -C(3,5);
    
    
else
    C = zeros(4,8);
    C(1,1) = -tan(gamma); C(1,2) = 1;
    C(2,1) = -tan(gamma); C(2,2) = -1;
    C(3,4) = -1;
    C(4,4) = 1;
end

