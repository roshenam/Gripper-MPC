function [loop] = Inv_Set_Calc(params)

omega = params.omega;
rp = params.rp; rs = params.rs;
Ts = params.Ts;
betaHIGH = params.betaHIGH;
eta = params.eta; tantol = params.tantol;
gamma = params.gamma;

% Creating dynamic system
A = zeros(6,6); A(1,4) = 1; A(2,5) = 1; A(3,6) = 1;
A(4,1) = omega^2; A(4,5) = 2*omega;
A(5,2) = omega^2; A(5,4) = -2*omega;
B = zeros(6,3); B(4,1) = 1; B(5,2) = 1; B(6,3) = 1;
C = zeros(6,6);  D = zeros(6,3);
sysD = ss(A,B,C,D);
sysD = c2d(sysD,Ts);

Ad = sysD.a; Bd = sysD.b; % Matrices for discrete system

% Y is a vector of constrained outputs. C and D are matrices defining Y's
% dependence on the states and control inputs

Dd = zeros(5,3);
%Dd = zeros(3,3);
Cd = zeros(5,6);
%Cd = zeros(2,6);
Cd(1,1) = -tan(gamma); Cd(1,2) = 1;
Cd(2,1) = -tan(gamma); Cd(2,2) = -1;
Cd(3,1) = -params.eta; Cd(3,2) = -params.eta; Cd(3,4) = -1;
Cd(4,1) = -params.eta; Cd(4,2) = -params.eta; Cd(4,5) = 1;
Cd(5,1) = -params.eta; Cd(5,2) = -params.eta; Cd(5,5) = -1;

% Upper bounds on thrust inputs and slack variables
Umax = params.Umax;
Tmax = params.Tmax; % max torque applied by flywheel
Ymax = [0; 0; -betaHIGH-eta*(rp+rs); tantol-eta*(rp+rs); tantol-eta*(rp+rs)];

% Creating weighting matrices for cost function
Q = params.Qval.*eye(6); Q(4:6,4:6) = .5*params.Qval.*eye(3); R = params.Rval.*eye(3);
system = LTISystem('A',Ad,'B',Bd,'C',Cd,'D',Dd,'Ts',Ts);
system.u.min = [-Umax; -Umax; -Tmax];
system.u.max = [Umax; Umax; Tmax];
system.y.max = Ymax;
system.x.penalty = QuadFunction(Q);
system.u.penalty = QuadFunction(R);
Pn = system.LQRPenalty;
Tset = system.LQRSet;
system.x.with('terminalPenalty');
system.x.terminalPenalty = Pn;
system.x.with('terminalSet');
system.x.terminalSet = Tset;

% Simulation
N = params.N;
mpc = MPCController(system,N);
loop = ClosedLoop(mpc, system);