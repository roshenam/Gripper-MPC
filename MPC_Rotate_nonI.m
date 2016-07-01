function [xtot, utot, cost, time, ytot] = MPC_Rotate_nonI(init,params,phi,omega)
% 2D accounting for platforms that rotate with a constant speed.
% Implementing cone constraints as soft through use of slack variables.
% This version is modeled in the non-inertial reference frame to make
% certain constraints constant

% Extracting initial conditions and parameters from inputs
x0 = init(1); y0 = init(2); %theta0 = init(3); 
vx0 = init(4); vy0 = init(5); %thetadot0 = init(6);
rp = params.rp; rs = params.rs; 
Ts = params.Ts; 
N = params.N; Nc = params.Nc; 
gamma = params.gamma;
Rmat = [cos(phi) sin(phi); -sin(phi) cos(phi)];
r0 = Rmat*[x0;y0];
v0 = Rmat*[vx0;vy0];
x0vec = [r0(1); r0(2); v0(1); v0(2); phi; omega];

% Creating dynamic system
A = zeros(6,6); A(1,3) = 1; A(2,4) = 1;
A(3,1) = omega^2; A(3,4) = 2*omega; 
A(4,2) = omega^2; A(4,3) = -2*omega;
A(5,6) = 1; 
B = zeros(6,4); B(3,1) = 1; B(4,2) = 1;
C = zeros(6,6);  D = zeros(6,4);
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

Dbig = zeros(2,4); Dbig(1,3) = 1; Dbig(2,4) = 1; % Slack variables to make constraints soft
Cbig = zeros(2,6);
Cbig(1,1) = -tan(gamma); Cbig(1,2) = 1;
Cbig(2,1) = -tan(gamma); Cbig(2,2) = -1;


% Upper bounds on thrust inputs and slack variables
Umax = params.Umax;
%Tmax = params.Tmax; % max torque applied by flywheel
Ymax = [0;0];

% Creating weighting matrices for cost function
Q = params.Qval.*eye(6); Q(5,5) = 0; Q(6,6) = 0; R = params.Rval.*eye(4); 
R(3,3) = params.slackweight; R(4,4) = params.slackweight; 
% Solving infinite horizon unconstrained LQR for stability enforcement
[~,P,~] = lqrd(A(1:4,1:4),B(1:4,1:2),Q(1:4,1:4),R(1:2,1:2),Ts);
Pbig = zeros(6,6);
Pbig(1:4,1:4) = P;

n = 6; m = 4; p=2; 
time = [];
utot = [];
ytot = [Ymax];
xtot = x0vec;
cost = [];
counter = 0;
%while norm(x0vec(1:2))>=(rp+rs)
for j=1:50
    counter = counter + 1;
    disp(['Running optimization ',num2str(counter),', distance from origin is ',num2str(norm(x0vec(1:2)))])
    tic
    cvx_begin quiet
    variables X(n,N+1) U(m,N) Y(p,N)
    X(:,2:N+1) == Abig*X(:,1:N) + Bbig*U;
    Y(:,1:N) == Cbig*X(:,1:N) + Dbig*U;
    X(:,1) == x0vec;
    max(Y(:,1:Nc)') <= Ymax';
    max((U(1,:).^2 + U(2,:).^2)') <= Umax';
    minimize (norm(Q*X(:,1:N),'fro') + norm(R*U(:,1:N),'fro') +...
        X(:,N+1)'*Pbig*X(:,N+1));
    cvx_end
    timecurr = toc;
    u = U(:,1);
    if isnan(u(1))
        disp(['Problem became infeasible at iteration ',num2str(counter)])
        break
    end
    x0vec = Abig*x0vec+Bbig*u;
    time = [time timecurr];
    xtot = [xtot x0vec] ;
    utot = [utot u];
    ytot = [ytot Ymax];
    currcost = cvx_optval;
    cost = [cost; currcost];
   
end
xs = xtot(1,:).*cos(xtot(5,:))-xtot(2,:).*sin(xtot(5,:));
ys = xtot(1,:).*sin(xtot(5,:))+xtot(2,:).*cos(xtot(5,:));
vxs = xtot(3,:).*cos(xtot(5,:))-xtot(4,:).*sin(xtot(5,:));
vys = xtot(3,:).*sin(xtot(5,:))+xtot(4,:).*cos(xtot(5,:));
xtot(7,:) = xtot(1,:); xtot(8,:) = xtot(2,:);
xtot(1,:) = xs; xtot(2,:) = ys;
xtot(9,:) = xtot(3,:); xtot(10,:) = xtot(4,:);
xtot(3,:) = vxs; xtot(4,:) = vys;
disp('Simulation Complete')
