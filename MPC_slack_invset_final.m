% Defining parameters
clear all
Matrix_Calc

phi = pi/4; omega = 2*pi/180;
params.phi = phi; % Phi is the angle of the port from the defined horizontal
params.omega = omega; % Omega is the rotation in rad/sec of the target
params.rp = .15; params.rtol = rtol;  params.rs = .15; params.gamma = gamma;
params.w = .01;
params.l = .1;
% rp is the radius of the target, rs is the radius of the spacecraft, gamma
% is the half angle of the constraint cone
params.Umax = h_u(1:5);
% Max torque in N*m/(kg*m^2)
params.Ts = T_MPC;
% Discretization constant (in seconds)
% Weights on states, control inputs, and slack variables
params.betaHIGH = -1.5; params.betaLOW = -0.2;
params.eta = 1;
params.tantol = .4;
% Constants defining the shape of the bounding function for normal
% component of hte velocity

rp = params.rp; rs = params.rs;
Ts = T_MPC;
betaHIGH = params.betaHIGH; betaLOW = params.betaLOW;
eta = params.eta; tantol = params.tantol;

Aw = zeros(4,4); Bw = zeros(4, 2);
Aw(1,3) = omega^2; Aw(1,2) = 2*omega;
Aw(2,4) = omega^2; Aw(2,1) = -2*omega;
C = zeros(4,4);  D = zeros(4,2);
sysD = ss(Aw,Bw,C,D);
sysD = c2d(sysD,Ts);

Awd = sysD.a;
Awd = [Awd zeros(4,2); zeros(2,6)]

Anew = A + Awd;

% Y is a vector of constrained outputs. C and D are matrices defining Y's
% dependence on the states and control inputs

Dd = zeros(2,5);
%Dd = zeros(3,3);
Cd = zeros(2,6);
%Cd = zeros(2,6);
Cd(1,1) = -tan(gamma); Cd(1,2) = 1; Dd(1, 4) = 1;
Cd(2,1) = -tan(gamma); Cd(2,2) = -1; Dd(2, 5) = 1;
%Cd(3,1) = -params.eta; Cd(3,2) = -params.eta; Cd(3,4) = -1;
%Cd(4,1) = -params.eta; Cd(4,2) = -params.eta; Cd(4,5) = 1;
%Cd(5,1) = -params.eta; Cd(5,2) = -params.eta; Cd(5,5) = -1;

%Cd(3,1) = eta1; Cd(3,2) = eta1; Cd(3,4) = -cos(gamma); Cd(3,5) = -sin(gamma);
%Cd(4,1) = -eta2; Cd(4,2) = -eta2; Cd(4,4) = cos(gamma); Cd(4,5) = sin(gamma);

K = -dlqr(Anew,B,Q,dlqr_controlweight*R);
eig(Anew+B*K)
print('second time')
P = dare(A+B*K, B, Q, dlqr_controlweight*R);
Q_tilde_f = P;          % p.d. weight matrix for final state

% Upper bounds on thrust inputs and slack variables
%Ymax = [0; 0; 0; 0];
%Ymax = [0; 0];
%Ymax = [tan(gamma)*params.rtol; tan(gamma)*params.rtol;...
    %-betaHIGH-eta*(rp+rs); tantol-eta*(rp+rs); tantol-eta*(rp+rs)];
%Ymax = [tan(gamma)*params.rtol; tan(gamma)*params.rtol];
%Ymax = [0; 0; -betaHIGH-eta*(rp+rs)];

model = LTISystem('A',Anew,'B',B(:,1:3));
model.x.min = [x1dot_min x2dot_min x1_min x2_min thetadot_min theta_min]';
model.x.max = [x1dot_max x2dot_max x1_max x2_max thetadot_max theta_max]';
model.u.min = [u1_min u2_min u3_min]';
model.u.max = [u1_max u2_max u3_max]';

%model.y.max = Ymax;
%system.y.with('softMax');
%TSet = model.invariantSet('maxIterations',40);
model.x.penalty = QuadFunction(Q_tilde);
model.u.penalty = QuadFunction(R(1:3,1:3));
Pn = model.LQRPenalty;
TSet = model.LQRSet
model.x.with('terminalPenalty');
model.x.terminalPenalty = Pn;
model.x.with('terminalSet');
model.x.terminalSet = TSet;
%%
clc
tic
params.N = 20;
N = params.N;
Nsim = 40;
params.Nsim = Nsim;
mpc = MPCController(model,N);
%exp = mpc.toExplicit();
loop = ClosedLoop(mpc, model);

% Specify initial conditions in frame with x axis intersecting with target
% point. Then rotate into shifted frame with x axis along bottom edge of
% cone
x0 =  3; y0 = 1; theta0 = 0; vx0 = 1; vy0 = 1; thetadot0 = 5*pi/180;
Rmat = [cos(phi) -sin(phi); sin(phi) cos(phi)];
rI = Rmat*[x0;y0];
vI = Rmat*[vx0;vy0];
init = [rI(1) rI(2) theta0 vI(1) vI(2) thetadot0];

%r0 = Rmat*[x0;y0];
%v0 = Rmat*[vx0;vy0];
nu0 = theta0-phi; nu0dot = thetadot0-omega;
%x0vec = [r0(1); r0(2); nu0; v0(1); v0(2); nu0dot];
x0vec = [vx0; vy0; x0; y0; nu0; nu0dot];
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
xtot(7,:) = phi+Ts*omega.*mult;
xtot(1,:) = xtotnonI(1,:).*cos(xtot(7,:))-xtotnonI(2,:).*sin(xtot(7,:));
xtot(2,:) = xtotnonI(1,:).*sin(xtot(7,:))+xtotnonI(2,:).*cos(xtot(7,:));
xtot(3,:) = xtotnonI(3,:) + xtot(7,:);
xtot(4,:) = xtotnonI(4,:).*cos(xtot(7,:))-xtotnonI(5,:).*sin(xtot(7,:));
xtot(5,:) = xtotnonI(4,:).*sin(xtot(7,:))+xtotnonI(5,:).*cos(xtot(7,:));
xtot(6,:) = xtotnonI(6,:) + omega;

utot(1,:) = utotnonI(1,:).*cos(xtot(7,1:end-1))-utotnonI(2,:).*sin(xtot(7,1:end-1));
utot(2,:) = utotnonI(1,:).*sin(xtot(7,1:end-1))+utotnonI(2,:).*cos(xtot(7,1:end-1));
utot(3,:) = utotnonI(3,:);
[vnorm, angattack, offset,collI] = Collision_compute(params,data)
Animate(init,params,xtot(:,1:collI),xtot(end,1:collI),0,'2D_Inv_Set',1,1)
% for i=1:length(data.X(1,:))
%     Collision_make(params, data.X(1:2,i), data.X(3,i), 1);
% end

%% Determining incoming conditions upon impact

% d_grip_offset = 0.156;                              % [m] Distance between gripper CG and target CG when captured
% d_grippad_offset = 0.041;                           % [m] Distance between gripper CG and gripper pad
% l_finger = 0.13;                                    % Length of gripper finger

