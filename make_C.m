function C = make_C(params, xk, yk, phik)

rp = params.rp; rtol = params.rtol; gamma = params.gamma;
C = zeros(2,6);
L1 = sin(phik+gamma)/((rp-rtol)*sin(gamma));
L2 = cos(phik+gamma)/((rp-rtol)*sin(gamma));
L3 = sin(phik-gamma)/((rp-rtol)*sin(gamma));
L4 = cos(phik-gamma)/((rp-rtol)*sin(gamma));

% Different than the paper. This is what I got from my math. Check later
% with their math
% C(1,1) = L1 - phik*L2;
% C(1,2) = -(L2 + phik*L1);
% C(1,10) = L1*yk + L2*xk;
% C(2,1) = -L3 - phik*L2;
% C(2,2) = L4 - phik*L3;
% C(2,10) = -L3*yk - L4*xk;

C(1,1) = L1;
C(1,2) = -L2;
C(2,1) = -L3;
C(2,2) = L4;
C(1,6) = L2*xk+L1*yk;
C(2,6) = -L4*xk-L3*yk;

