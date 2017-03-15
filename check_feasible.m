 
Cd(1,1) = -tan(pi/4+gamma); Cd(1,2) = 1; Dd(1,4) = 1; % LOS constraints
Cd(2,1) = tan(pi/4-gamma); Cd(2,2) = -1; Dd(2,5) = 1;
 Cd(3,1) = -eta_vel; Cd(3,2) = -eta_vel; Cd(3,4) = -1; % Vn upper bound
 Cd(4,1) = -eta_vel; Cd(4,2) = -eta_vel; Cd(4,4) = 1;
 Cd(5,1) = -eta_vel; Cd(5,2) = -eta_vel; Cd(5,5) = 1; % Vn constraints
 Cd(6,1) = -eta_vel; Cd(6,2) = -eta_vel; Cd(6,5) = -1;
% %Cd(6,1) = -eta_theta; Cd(6,2) = -eta_theta; Cd(6,3) = 1; % Theta constraints
% %Cd(7,1) = -eta_theta; Cd(7,2) = -eta_theta; Cd(7,3) = -1;
% %Cd(6,1) = -1; Cd(6,2) = -1; 
 Cd(7,1) = -eta_omega; Cd(7,2) = -eta_omega; Cd(7,6) = 1; % Omega constraints 
 Cd(8,1) = -eta_omega; Cd(8,2) = -eta_omega; Cd(8,6) = -1;
 Cd(9,4) = 1; % Restricting to 3rd quadrant
 Cd(10,5) = 1; 
 Cd(11,1) = -1;
 Cd(12,1) = -1;

%Cd(3,1) = eta1; Cd(3,2) = eta1; Cd(3,4) = -cos(gamma); Cd(3,5) = -sin(gamma);
%Cd(4,1) = -eta2; Cd(4,2) = -eta2; Cd(4,4) = cos(gamma); Cd(4,5) = sin(gamma);

% Upper bounds on thrust inputs and slack variables
Umax = params.Umax;
Tmax = params.Tmax; % max torque applied by flywheel
'hi'
%Ymax = [0; 0; 0; 0];
%Ymax = [0; 0];
%Ymax = [tan(gamma)*params.rtol; tan(gamma)*params.rtol;...
   % -betaHIGH-eta*(rp+rs); tantol-eta*(rp+rs); tantol-eta*(rp+rs)];
%Ymax = -betaHIGH-eta*(rp+rs);
%Ymax = [0; 0; -betaHIGH-eta*(rp+rs)];
Ymax = [(params.rp-params.rtol)/(1-tan(pi/4+gamma)); -(params.rp-params.rtol)/(1-tan(pi/4-gamma));...]%; tan(pi/4-gamma)*(params.rp-params.rtol)]%;...
    0; 0; 0; 0; 0; 0; 0; 0; 0; 0]%; 0]%; 100*vmax/sqrt(2)];...
     %vmax/sqrt(2)-eta*(rp+rs)]%; vmax/sqrt(2)-eta*(rp+rs); vmax/sqrt(2)-eta*(rp+rs);...
%     theta_c - eta*(rp+rs); theta_c - eta*(rp+rs); -eta*(rp+rs);...
%     omega_max - eta*(rp+rs); omega_max - eta*(rp+rs)];
   % -betaHIGH-eta*(rp+rs); tantol-eta*(rp+rs); tantol-eta*(rp+rs)];

% Creating weighting matrices for cost function
Q = params.Qval.*eye(6); Q(3,3) = 10000;
%R = params.Rval.*eye(3); R(3,3) = params.Rval/2;
R = params.Rval.*eye(3); R(3,3) = 10; R(4,4) = 10^6; R(5,5) = 10^6;
system = LTISystem('A',Ad,'B',Bd,'C',Cd,'D',Dd,'Ts',Ts);
% system.u.min = [-Umax; -Umax; -Tmax];
% system.u.max = [Umax; Umax; Tmax];

system.u.min = [-Umax; -Umax; -Tmax; -10^4; -10^4];
system.u.max = [Umax; Umax; Tmax; 10^4; 10^4];
system.y.max = Ymax;
%system.x.min = [0; -inf; -inf; -inf; -inf; -inf];
%system.y.with('softMax');
%xref = [0 0 0 -.2 0 0 0];
system.x.penalty = QuadFunction(Q);
system.u.penalty = QuadFunction(R);
Pn = system.LQRPenalty;
Tset = system.LQRSet;
%Tset = system.invariantSet('maxiterations',40)
system.x.with('terminalPenalty');
system.x.terminalPenalty = Pn;
system.x.with('terminalSet');
system.x.terminalSet = Tset;

clc
tic
params.N = 30;
N = params.N;
Nsim = 40;
params.Nsim = Nsim;
mpc = MPCController(system,N);
%exp = mpc.toExplicit();
loop = ClosedLoop(mpc, system);

% Specify initial conditions in frame with x axis intersecting with target
% point. Then rotate into shifted frame with x axis along bottom edge of
% cone
x0 =  1; y0 = 3; theta0 = 0; vx0 = -.05; vy0 = -.05; thetadot0 = 0;
Rmat = [cos(phi) -sin(phi); sin(phi) cos(phi)];
rI = Rmat*[x0;y0];
vI = Rmat*[vx0;vy0];
init = [rI(1) rI(2) theta0 vI(1) vI(2) thetadot0];

%r0 = Rmat*[x0;y0];
%v0 = Rmat*[vx0;vy0];
nu0 = theta0-phi+pi/4; nu0dot = thetadot0-omega;
%x0vec = [r0(1); r0(2); nu0; v0(1); v0(2); nu0dot];
x0vec = [x0; y0; nu0; vx0; vy0; nu0dot];

try
    data = loop.simulate(x0vec(1:6), Nsim);
catch
    check = -1
    error_flag = 1;
    'hi'
end

if error_flag
else
    if size(data.X,2) == 1
        check = -1;
        'bye'
    else
        check = 1;
    end
end