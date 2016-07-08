% Defining parameters 
phi = pi/4; omega = 5*pi/180;
params.phi = phi;
params.omega = omega;
params.rp = .2; params.rtol = .2+.01;  params.rs = .2; params.gamma = pi/20;
params.Umax = 1;
params.Tmax = 1;
params.Ts = .2;
params.Nc = 15;
params.Qval = 10^5; params.Rval = 10^2; params.slackweight = 10^8;
params.betaHIGH = 1.5; params.betaLOW = 0.2;

rp = params.rp; rs = params.rs; 
Ts = params.Ts; 
eta = params.eta; betaHIGH = params.betaHIGH; betaLOW = params.betaLOW;
gamma = params.gamma;
eta1 = betaHIGH/((rp+rs)*(cos(gamma)+sin(gamma))); 
eta2 = betaLOW/((rp+rs)*(cos(gamma)+sin(gamma))); 

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

Dd = zeros(4,3); 
Cd = zeros(4,6);
Cd(1,1) = -tan(gamma); Cd(1,2) = 1;
Cd(2,1) = -tan(gamma); Cd(2,2) = -1;
Cd(3,1) = -eta1; Cd(3,2) = -eta1; Cd(3,4) = -cos(gamma); Cd(3,5) = -sin(gamma);
Cd(4,1) = -eta2; Cd(4,2) = -eta2; Cd(4,4) = cos(gamma); Cd(4,5) = sin(gamma);



% Upper bounds on thrust inputs and slack variables
Umax = params.Umax;
Tmax = params.Tmax; % max torque applied by flywheel
Ymax = [0; 0; 0; 0];

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
mpc = MPCController(system,N);
%exp = mpc.toExplicit();
loop = ClosedLoop(mpc, system);
params.N = 30;
N = params.N; 
Nsim = 17;
params.Nsim = Nsim;

x0 =  5; y0 = 5; theta0 = -pi/4; vx0 = -1.5; vy0 = -1.5; thetadot0 = 50*pi/180;
init = [x0 y0 theta0 vx0 vy0 thetadot0];
Rmat = [cos(phi-gamma) sin(phi-gamma); -sin(phi-gamma) cos(phi-gamma)];
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
xtotnonI = data.X; % Save data in non-inertial frame
utotnonI = data.U; % Save thrust data in non-inertial frame
xtot = zeros(8,Nsim+1);
utot = zeros(3,Nsim);
mult = 0:Nsim;
xtot(8,:) = phi+Ts*omega.*mult;
xtot(1,:) = xtotnonI(1,:).*cos(xtot(8,:)-gamma)-xtotnonI(2,:).*sin(xtot(8,:)-gamma);
xtot(2,:) = xtotnonI(1,:).*sin(xtot(8,:)-gamma)+xtotnonI(2,:).*cos(xtot(8,:)-gamma);
xtot(4,:) = xtotnonI(4,:).*cos(xtot(8,:)-gamma)-xtotnonI(5,:).*sin(xtot(8,:)-gamma);
xtot(5,:) = xtotnonI(4,:).*sin(xtot(8,:)-gamma)+xtotnonI(5,:).*cos(xtot(8,:)-gamma);
xtot(3,:) = xtotnonI(3,:) + xtot(8,:);
xtot(6,:) = xtotnonI(6,:) + omega;
xtot(9,:) = eta1.*xtotnonI(1,:)+eta1.*xtotnonI(2,:);
xtot(10,:) = eta2.*xtotnonI(1,:)+eta2.*xtotnonI(2,:);
xtot(11,:) = xtotnonI(4,:).*cos(gamma)+xtotnonI(5,:).*sin(gamma);
utot(1,:) = utotnonI(1,:).*cos(xtot(8,1:end-1)-gamma)-utotnonI(2,:).*sin(xtot(8,1:end-1)-gamma);
utot(2,:) = utotnonI(1,:).*sin(xtot(8,1:end-1)-gamma)+utotnonI(2,:).*cos(xtot(8,1:end-1)-gamma);
utot(3,:) = utotnonI(3,:);
%InvSet = system.invariantSet();
%pause