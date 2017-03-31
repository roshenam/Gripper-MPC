% Defining parameters
clear all
phi = pi/4; omega = 2*pi/180;
params.phi = phi; % Phi is the angle of the port from the defined horizontal
params.omega = omega; % Omega is the rotation in rad/sec of the target
params.rp = .15; params.rtol = params.rp+.02;  params.rs = .15; params.gamma = 50*pi/180;
%params.rtol = params.rp;
params.w = .01;
params.l = .1;
%theta_c = 4*pi/180;
theta_c = 2*pi/180;
gamma_max = atan(.05/(params.rp+params.rs)) - theta_c
%gamma_max = 5*pi/180;
params.gamma = gamma_max;
params.mass = 18.08; %[kg]
%%
omega_c = 100*pi/180;

% rp is the radius of the target, rs is the radius of the spacecraft, gamma
% is the half angle of the constraint cone
%params.Umax = .2;
params.Umax = .4/params.mass;
% Max thrust in Newtons/kg of the free flyer
params.Tmax = .3; % 50% of stall torque
% Max torque in N*m/(kg*m^2)
params.Ts = 0.5;
% Discretization constant (in seconds)
params.Qval = 10^3; params.Rval = 10^2;
params.Qval2 = 10^2;
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
vmax = .3; % velocity used on d and omega plots from old paper
% Creating dynamic system
A = zeros(6,6); A(1,4) = 1; A(2,5) = 1; A(3,6) = 1;
A(4,1) = omega^2; A(4,5) = 2*omega;
A(5,2) = omega^2; A(5,4) = -2*omega;
B = zeros(6,5); B(4,1) = 1; B(5,2) = 1; B(6,3) = 1;
C = zeros(6,6);  D = zeros(6,5);
sysD = ss(A,B,C,D);
sysD = c2d(sysD,Ts);

Ad = sysD.a; Bd = sysD.b; % Matrices for discrete system

% Y is a vector of constrained outputs. C and D are matrices defining Y's
% dependence on the states and control inputs

%Dd = zeros(5,3);
%Dd = zeros(3,3);
Dd = zeros(10,5);
Cd = zeros(10,6);
%Cd = zeros(2,6);
eta_vel = vmax/(sqrt(2)*(rp+rs));
eta_theta = theta_c/(rp+rs);
eta_omega = omega_c/(rp+rs);
Cd(1,1) = -tan(pi/4+gamma); Cd(1,2) = 1; Dd(1,4) = 1; % LOS constraints
Cd(2,1) = tan(pi/4-gamma); Cd(2,2) = -1; Dd(2,5) = 1;
 Cd(3,1) = -eta_vel; Cd(3,2) = -eta_vel; Cd(3,4) = -1; % Vn upper bound
 Cd(4,1) = -eta_vel; Cd(4,2) = -eta_vel; Cd(4,4) = 1;
 Cd(5,1) = -eta_vel; Cd(5,2) = -eta_vel; Cd(5,5) = 1; % Vn constraints
 Cd(6,1) = -eta_vel; Cd(6,2) = -eta_vel; Cd(6,5) = -1;
Cd(7,4) = 1; % Angle of attack constraints through velocity 
Cd(8,5) = 1; 
Cd(9,1) = -1; % x > 0
Cd(10,2) = -1; % y > 0 

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
    0; 0; 0; 0; 0; 0; 0; 0]%; 0]%; 100*vmax/sqrt(2)];...
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
%system.x.with('reference')
%xref = [0 0 0 -.2 0 0 0];
system.x.penalty = QuadFunction(Q);
system.u.penalty = QuadFunction(R);
Pn = system.LQRPenalty;
%Tset = system.LQRSet
Tset = system.invariantSet()
%Tset = system.invariantSet('maxiterations',40)
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

% Specify initial conditions in frame with target point at pi/4
x0 =  1; y0 = 1.5; theta0 = phi+pi/4; vx0 = 0; vy0 = 0; thetadot0 = 0;
Rmat = [cos(phi) -sin(phi); sin(phi) cos(phi)];
rI = Rmat*[x0;y0];
vI = Rmat*[vx0;vy0];
init = [rI(1) rI(2) theta0 vI(1) vI(2) thetadot0];

%r0 = Rmat*[x0;y0];
%v0 = Rmat*[vx0;vy0];
nu0 = theta0-phi-pi/4; nu0dot = thetadot0-omega;
%x0vec = [r0(1); r0(2); nu0; v0(1); v0(2); nu0dot];
x0vec = [x0; y0; nu0; vx0; vy0; nu0dot];
%data = loop.simulate(x0vec(1:6), Nsim,'x.reference',xref')
data = loop.simulate(x0vec(1:6), Nsim)

toc
xtotnonI = data.X; % Save data in non-inertial frame
utotnonI = data.U; % Save thrust data in non-inertial frame
vnbound = -eta.*data.X(1,:)-eta.*data.X(2,:)+betaHIGH+eta*(rp+rs);
vtbound = -eta.*data.X(1,:)-eta.*data.X(2,:)-tantol+eta*(rp+rs);
xtot = zeros(7,Nsim+1);
utot = zeros(3,Nsim);
mult = 0:Nsim;
time = Ts.*mult;
xtot(7,:) = phi+Ts*omega.*mult-pi/4;
xtot(8,:) = xtot(7,:)+pi/4;
xtot(1,:) = xtotnonI(1,:).*cos(xtot(7,:))-xtotnonI(2,:).*sin(xtot(7,:));
xtot(2,:) = xtotnonI(1,:).*sin(xtot(7,:))+xtotnonI(2,:).*cos(xtot(7,:));
xtot(3,:) = xtotnonI(3,:) + xtot(8,:);
xtot(4,:) = xtotnonI(4,:).*cos(xtot(7,:))-xtotnonI(5,:).*sin(xtot(7,:));
xtot(5,:) = xtotnonI(4,:).*sin(xtot(7,:))+xtotnonI(5,:).*cos(xtot(7,:));
xtot(6,:) = xtotnonI(6,:) + omega;

utot(1,:) = utotnonI(1,:).*cos(xtot(7,1:end-1))-utotnonI(2,:).*sin(xtot(7,1:end-1));
utot(2,:) = utotnonI(1,:).*sin(xtot(7,1:end-1))+utotnonI(2,:).*cos(xtot(7,1:end-1));
utot(3,:) = utotnonI(3,:);
[vnorm, angattack, offset,collI] = Collision_compute(params,data)
%collI = length(data.X);
%Animate(init,params,xtot(:,1:collI),xtot(end,1:collI),0,'2D_Inv_Set',1,1)
Animate(init,params,xtot(:,1:collI-1),xtot(end,1:collI-1),0,'2D_Inv_Set',1,1)
%Animate(init,params,xtot(:,:),xtot(end,:),0,'2D_Inv_Set',1,1)
% for i=1:length(data.X(1,:))
%     Collision_make(params, data.X(1:2,i), data.X(3,i), 1);
% end

%% Determining incoming conditions upon impact

% d_grip_offset = 0.156;                              % [m] Distance between gripper CG and target CG when captured
% d_grippad_offset = 0.041;                           % [m] Distance between gripper CG and gripper pad
% l_finger = 0.13;                                    % Length of gripper finger

