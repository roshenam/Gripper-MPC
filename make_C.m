function C = make_C(params, x0vec, inertial)

rp = params.rp; rtol = params.rtol; gamma = params.gamma;
xk = x0vec(1); yk = x0vec(2);
phik = x0vec(8);
if inertial == 1
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
    C(1,8) = L2*xk+L1*yk;
    C(2,8) = -L4*xk-L3*yk;
    
    % Normal velocity constraints
    C(3,4) = -cos(phik);
    C(3,5) = -sin(phik);
    C(4,4) = -C(3,4);
    C(4,5) = -C(3,5);
    
    
elseif inertial == 0
    % Upper and lower bounds on normal velocity 
    C = zeros(4,8);
    C(1,1) = -tan(gamma); C(1,2) = 1;
    C(2,1) = -tan(gamma); C(2,2) = -1;
    C(3,4) = -1;
    C(4,4) = 1;
    
elseif inertial == 2
    %C = zeros(4,10);
    C = zeros(2,10);
    C(1,1) = -tan(gamma); C(1,2) = 1; C(1,7) = -tan(gamma)*params.vD;
    C(2,1) = -tan(gamma); C(2,2) = -1; C(2,7) = -tan(gamma)*params.vD;
    %C(3,1) = -params.eta; C(3,2) = -params.eta; 
    %C(3,4) = 1; C(3,7) = -params.eta*params.vD;
    %C(4,1) = -params.eta; C(4,2) = -params.eta; 
    %C(4,4) = -1; C(4,7) = -params.eta*params.vD; 
elseif inertial == 3
    % Upper bound on normal velocity, upper and lower bounds on tangential
    % velocity
    C = zeros(7,8);
    C(1,1) = -tan(gamma); C(1,2) = 1;
    C(2,1) = -tan(gamma); C(2,2) = -1;
    C(3,4) = -1; % upper bound on normal velocity
    C(4,5) = 1; C(5,5) = -1; % upper bound on tangential velocity
    C(6,3) = 1; C(7,3) = -1; % bounds on omega
end

