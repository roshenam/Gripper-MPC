% Defining parameters
phi = pi/4; omega = 2*pi/180;
params.phi = phi; % Phi is the angle of the port from the defined horizontal
params.omega = omega; % Omega is the rotation in rad/sec of the target
params.rp = .2; params.rtol = .2+.01;  params.rs = .2; params.gamma = 20*pi/180;
% rp is the radius of the target, rs is the radius of the spacecraft, gamma
% is the half angle of the constraint cone
params.Umax = .5;
% Max thrust in Newtons/kg of the free flyer
params.Tmax = 1;
% Max torque in N*m/(kg*m^2)
params.Ts = .2;
% Discretization constant (0.2 seconds)
params.Nc = 15;
% Control Horizon for MPC
params.Qval = 10^4; params.Rval = 10^4;
% Weights on states, control inputs, and slack variables
params.betaHIGH = -1.5; params.betaLOW = -0.2;
params.eta = 1;
params.tantol = .4;
% Constants defining the shape of the bounding function for normal
% component of hte velocity

rp = params.rp; rs = params.rs;
Ts = params.Ts;
betaHIGH = params.betaHIGH; betaLOW = params.betaLOW;
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

%Cd(3,1) = eta1; Cd(3,2) = eta1; Cd(3,4) = -cos(gamma); Cd(3,5) = -sin(gamma);
%Cd(4,1) = -eta2; Cd(4,2) = -eta2; Cd(4,4) = cos(gamma); Cd(4,5) = sin(gamma);



% Upper bounds on thrust inputs and slack variables
Umax = params.Umax;
Tmax = params.Tmax; % max torque applied by flywheel
%Ymax = [0; 0; 0; 0];
%Ymax = [0; 0];
Ymax = [0; 0; -betaHIGH-eta*(rp+rs); tantol-eta*(rp+rs); tantol-eta*(rp+rs)];
%Ymax = [0; 0; -betaHIGH-eta*(rp+rs)];

% Creating weighting matrices for cost function
Q = params.Qval.*eye(6); Q(4:6,4:6) = .5*params.Qval.*eye(3); R = params.Rval.*eye(3);
system = LTISystem('A',Ad,'B',Bd,'C',Cd,'D',Dd,'Ts',Ts);
system.u.min = [-Umax; -Umax; -Tmax];
system.u.max = [Umax; Umax; Tmax];
system.y.max = Ymax;
%system.y.with('softMax');
system.x.penalty = QuadFunction(Q);
system.u.penalty = QuadFunction(R);
Pn = system.LQRPenalty;
Tset = system.LQRSet;
system.x.with('terminalPenalty');
system.x.terminalPenalty = Pn;
system.x.with('terminalSet');
system.x.terminalSet = Tset;
%%
clc
tic
params.N = 30;
N = params.N;
Nsim = 30;
params.Nsim = Nsim;
mpc = MPCController(system,N);
%exp = mpc.toExplicit();
loop = ClosedLoop(mpc, system);

% Specify initial conditions in frame with x axis intersecting with target
% point. Then rotate into shifted frame with x axis along bottom edge of
% cone
x0 =  5; y0 = 0; theta0 = -pi/4; vx0 = -.9; vy0 = .9; thetadot0 = 5*pi/180;
Rmat = [cos(phi) -sin(phi); sin(phi) cos(phi)];
rI = Rmat*[x0;y0];
vI = Rmat*[vx0;vy0];
init = [rI(1) rI(2) theta0 vI(1) vI(2) thetadot0];

%r0 = Rmat*[x0;y0];
%v0 = Rmat*[vx0;vy0];
nu0 = theta0-phi; nu0dot = thetadot0-omega;
%x0vec = [r0(1); r0(2); nu0; v0(1); v0(2); nu0dot];
x0vec = [x0; y0; nu0; vx0; vy0; nu0dot];
data = loop.simulate(x0vec(1:6), Nsim);
toc
xtotnonI = data.X; % Save data in non-inertial frame
utotnonI = data.U; % Save thrust data in non-inertial frame
vnbound = -eta.*data.X(1,:)-eta.*data.X(2,:)+betaHIGH+eta*(rp+rs);
vtbound = -eta.*data.X(1,:)-eta.*data.X(2,:)-tantol+eta*(rp+rs);
xtot = zeros(8,Nsim+1);
utot = zeros(3,Nsim);
mult = 0:Nsim;
time = Ts.*mult;
xtot(8,:) = phi+Ts*omega.*mult;
xtot(1,:) = xtotnonI(1,:).*cos(xtot(8,:))-xtotnonI(2,:).*sin(xtot(8,:));
xtot(2,:) = xtotnonI(1,:).*sin(xtot(8,:))+xtotnonI(2,:).*cos(xtot(8,:));
xtot(3,:) = xtotnonI(3,:) + xtot(8,:);
xtot(4,:) = xtotnonI(4,:).*cos(xtot(8,:))-xtotnonI(5,:).*sin(xtot(8,:));
xtot(5,:) = xtotnonI(4,:).*sin(xtot(8,:))+xtotnonI(5,:).*cos(xtot(8,:));
xtot(6,:) = xtotnonI(6,:) + omega;

utot(1,:) = utotnonI(1,:).*cos(xtot(8,1:end-1))-utotnonI(2,:).*sin(xtot(8,1:end-1));
utot(2,:) = utotnonI(1,:).*sin(xtot(8,1:end-1))+utotnonI(2,:).*cos(xtot(8,1:end-1));
utot(3,:) = utotnonI(3,:);

Animate(init,params,xtot,xtot(end,:),0,'hi',1,1)
%% Determining incoming conditions upon impact

d_grip_offset = 0.156;                              % [m] Distance between gripper CG and target CG when captured
d_grippad_offset = 0.041;                           % [m] Distance between gripper CG and gripper pad
l_finger = 0.13;                                    % Length of gripper finger

idx = find(dist<=(rp+rs),1);
final.vn = data.X(4,idx-1);
final.vt = data.X(5,idx-1);
final.vmag = norm([data.X(4,idx-1) data.X(5,idx-1)]);
final.nu = data.X(3,idx-1)*180/pi;
final.angattack = atan2(data.X(5,idx-1),data.X(4,idx-1));
final.nudot = data.X(6,idx-1)*180/pi;
final.thrust = sum(sqrt(utot(1,:).^2+utot(2,:).^2));
d_rem = dist(idx-1) - (rp+rs);
t_rem = d_rem/final.vmag;

final