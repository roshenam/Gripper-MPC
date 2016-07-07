x0 =  3; y0 = 3; theta0 = -pi/4; vx0 = .2; vy0 = .2; thetadot0 = -20*pi/180;
init = [x0 y0 theta0 vx0 vy0 thetadot0];
phi = pi/4; omega = 5*pi/180;
params.phi = phi;
params.omega = omega;
params.Ro = 850*10^3; % [m] orbital radius of target
params.mu = 3.986004418*10^14; % [m^3/s^2] standard gravitational parameter of earth
params.rp = .2; params.rtol = .2+.01;  params.rs = .2; params.gamma = pi/20;
params.Umax = 2;
params.Tmax = 1;
params.Ts = .2;
params.N = 30; params.Nc = 15;
params.Qval = 10^5; params.Rval = 10^2; params.slackweight = 10^8;
params.eta = 1; params.betaHIGH = 1.5; params.betaLOW = 0.2;

x0 = init(1); y0 = init(2); theta0 = init(3); 
vx0 = init(4); vy0 = init(5); thetadot0 = init(6);
rp = params.rp; rs = params.rs; 
Ts = params.Ts; 
N = params.N; Nc = params.Nc; 
gamma = params.gamma;
Rmat = [cos(phi) sin(phi); -sin(phi) cos(phi)];
r0 = Rmat*[x0;y0];
v0 = Rmat*[vx0;vy0];
nu0 = theta0-phi; nu0dot = thetadot0-omega;
x0vec = [r0(1); r0(2); nu0; v0(1); v0(2); nu0dot; phi; omega];

% Creating dynamic system
A = zeros(8,8); A(1,4) = 1; A(2,5) = 1; A(3,6) = 1;
A(4,1) = omega^2; A(4,5) = 2*omega; 
A(5,2) = omega^2; A(5,4) = -2*omega;
A(7,8) = 1; 
B = zeros(8,5); B(4,1) = 1; B(5,2) = 1; B(6,3) = 1;
C = zeros(8,8);  D = zeros(8,5);
sysD = ss(A,B,C,D);
sysD = c2d(sysD,Ts);

% Adding additional states to the system to account for rotation of the
% platform. States will be rx and ry (the x and y pos of the selected
% point), and sigx and sigy (the x and y distance between the selected
% point and the center of mass of the spacecraft). Also adding states to
% predict the future state of the LOS cone constraints. State vector is now
% [x y x' y' z1 z2 theta theta']
% Abig = zeros(8,8); Abig(1:4,1:4) = sysD.a(1:4,1:4); 
% Abig(5,5) = 2; Abig(5,6) = -1; Abig(6,5) = 1;
% Abig(7:8,7:8) = sysD.a(5:6,5:6);
% Bbig = zeros(8,7); Bbig(1:4,1:5) = sysD.b(1:4,1:5);
% Bbig(7:8,1:5) = sysD.b(5:6,1:5);
Abig = sysD.a; Bbig = sysD.b;

% Y is a vector of constrained outputs. C and D are matrices defining Y's
% dependence on the states and control inputs

Dbig = zeros(2,5); Dbig(1,4) = 1; Dbig(2,5) = 1; % Slack variables to make constraints soft
Cbig = zeros(2,8);
Cbig(1,1) = -tan(gamma); Cbig(1,2) = 1;
Cbig(2,1) = -tan(gamma); Cbig(2,2) = -1;


% Upper bounds on thrust inputs and slack variables
Umax = params.Umax;
Tmax = params.Tmax; % max torque applied by flywheel
Ymax = [0;0];

% Creating weighting matrices for cost function
Q = params.Qval.*eye(8); Q(7,7) = 0; Q(8,8) = 0; R = params.Rval.*eye(5); 
R(4,4) = params.slackweight; R(5,5) = params.slackweight; 
% Solving infinite horizon unconstrained LQR for stability enforcement
[~,P,~] = lqrd(A(1:6,1:6),B(1:6,1:3),Q(1:6,1:6),R(1:3,1:3),Ts);
Pbig = zeros(8,8);
Pbig(1:6,1:6) = P;

n = 8; m = 5; p=2; 
time = [];
utot = [];
ytot = Ymax;
xtot = x0vec;
cost = [];
counter = 0;
system = LTISystem('A',Abig(1:6,1:6),'B',Bbig(1:6,1:3),'C',Cbig(:,1:6),'D',Dbig(:,1:3),'Ts',Ts);
system.u.min = [-Umax; -Umax; -Umax];
system.u.max = [Umax; Umax; Umax];
system.y.max = [0; 0];
system.x.penalty = QuadFunction(Q(1:6,1:6));
system.u.penalty = QuadFunction(R(1:3,1:3));
Pn = system.LQRPenalty;
Tset = system.LQRSet;
system.x.with('terminalPenalty');
system.x.terminalPenalty = Pn;
system.x.with('terminalSet');
system.x.terminalSet = Tset;
%%
clc
N = 10;
mpc = MPCController(system,N);
%exp = mpc.toExplicit();
loop = ClosedLoop(mpc, system);
Nsim = 20;
x0 =  1; y0 = 1; theta0 = -pi/4; vx0 = 1.5; vy0 = .5; thetadot0 = -100*pi/180;
r0 = Rmat*[x0;y0];
v0 = Rmat*[vx0;vy0];
nu0 = theta0-phi; nu0dot = thetadot0-omega;
x0vec = [r0(1); r0(2); nu0; v0(1); v0(2); nu0dot];
data = loop.simulate(x0vec(1:6), Nsim);
figure(1)
subplot(3,1,1)
plot(data.X(1,:),data.X(2,:),'k','linewidth',2)
xlabel('X'); ylabel('Y');
subplot(3,1,2)
plot(data.X(3,:),'k','linewidth',2)
xlabel('T'); ylabel('\theta');
subplot(3,1,3)
plot(data.U(1,:),'linewidth',2)
hold all
plot(data.U(2,:),'linewidth',2)
plot(data.U(3,:),'linewidth',2)
legend('Ux','Uy','Uz')
%InvSet = system.invariantSet();
%pause