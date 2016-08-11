function [xtot, utot, xtotnonI, utotnonI, final] = Inv_Set_Sim(x0vec, params)
% Defining parameters
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
Nsim = params.Nsim;
mpc = MPCController(system,N);
loop = ClosedLoop(mpc, system);
data = loop.simulate(x0vec(1:6), Nsim);
[~,n] = size(data.X);

if n==1
    xtot = NaN;
    utot = NaN;
    xtotnonI = NaN; 
    utotnonI =  NaN;
    final = NaN;
    return
end
    
xtotnonI = data.X; % Save data in non-inertial frame
utotnonI = data.U; % Save thrust data in non-inertial frame
xtot = zeros(7,Nsim+1);
utot = zeros(3,Nsim);
mult = 0:Nsim;

phi = params.phi;
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

% Determining incoming conditions upon impact
dist = sqrt(data.X(1,:).^2+data.X(2,:).^2);
idx = find(dist<=(rp+rs),1);
final.idx = idx;
final.vn = data.X(4,idx-1);
final.vt = data.X(5,idx-1);
final.vmag = norm([data.X(4,idx-1) data.X(5,idx-1)]);
final.nu = data.X(3,idx-1)*180/pi;
final.angattack = atan2(data.X(5,idx-1),data.X(4,idx-1));
final.nudot = data.X(6,idx-1)*180/pi;
final.thrust = sum(sqrt(utot(1,:).^2+utot(2,:).^2));

