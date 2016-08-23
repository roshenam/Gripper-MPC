function [xtot, xtotnonI, utot, cost, time, ytot, slack] = MPC_Rotate_nonI(init,params)
% 2D accounting for platforms that rotate with a constant speed.
% Implementing cone constraints as soft through use of slack variables.
% This version is modeled in the non-inertial reference frame to make
% constraints constant

% Extracting initial conditions and parameters from inputs
x0 = init(1); y0 = init(2); theta0 = init(3); 
vx0 = init(4); vy0 = init(5); thetadot0 = init(6);
rp = params.rp; rs = params.rs; 
Ts = params.Ts; 
N = params.N; Nc = params.Nc; 
phi = params.phi;
omega = params.omega;
Rmat = [cos(phi) sin(phi); -sin(phi) cos(phi)];
r0 = Rmat*[x0;y0];
v0 = Rmat*[vx0;vy0];
nu0 = theta0-phi; nu0dot = thetadot0-omega;
x0vec = [r0(1); r0(2); nu0; v0(1); v0(2); nu0dot; phi+omega*Ts; phi];

% Creating dynamic system
A = zeros(6,6); A(1,4) = 1; A(2,5) = 1; A(3,6) = 1;
A(4,1) = omega^2; A(4,5) = 2*omega; 
A(5,2) = omega^2; A(5,4) = -2*omega;
B = zeros(6,5); B(4,1) = 1; B(5,2) = 1; B(6,3) = 1;
C = zeros(6,6);  D = zeros(6,5);
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
Abig = zeros(8,8); Abig(1:6,1:6) = sysD.a; 
Abig(7,7) = 2; Abig(7,8) = -1; Abig(8,7) = 1;
Bbig = zeros(8,7); Bbig(1:6,1:3) = sysD.b(1:6,1:3);

% Y is a vector of constrained outputs. C and D are matrices defining Y's
% dependence on the states and control inputs

Dbig = zeros(2,7); Dbig(1,4) = 1; Dbig(2,5) = 1; % Slack variables to make constraints soft
Dbig(3,6) = 1; Dbig(4,7) = 1;
Cbig = make_C(params, x0vec, 0);

% Upper bounds on thrust inputs and slack variables
Umax = params.Umax;
Tmax = params.Tmax; % max torque applied by flywheel
Ymax = make_Ymax(params, x0vec, 0);

% Creating weighting matrices for cost function
Q = params.Qval.*eye(8); Q(7,7) = 0; Q(8,8) = 0; R = params.Rval.*eye(5);
if params.slack_variable
    dist = sqrt((r0(1)-rp*cos(nu0))^2+(r0(2)-rp*sin(nu0))^2);
    s_weight = 10^(params.slack_slope*dist + params.slack_intercept);
    R(4,4) = s_weight; R(5,5) = s_weight;
    R(6,6) = s_weight; R(7,7) = s_weight;
    slack = s_weight;
else
    
    R(4,4) = params.slackweight; R(5,5) = params.slackweight;
    R(6,6) = params.slackweight; R(7,7) = params.slackweight;
    slack = [];
end
% Solving infinite horizon unconstrained LQR for stability enforcement
[~,P,~] = lqrd(A(1:6,1:6),B(1:6,1:3),Q(1:6,1:6),R(1:3,1:3),Ts);
Pbig = zeros(8,8);
Pbig(1:6,1:6) = P;

n = 8; m = 7; p=4; 
time = [];
utot = [];
ytot = Ymax;
xtot = x0vec;
cost = [];
counter = 0;
while norm(x0vec(1:2))>=(rp+rs)
%for j=1:50
    counter = counter + 1;
    disp(['Running optimization ',num2str(counter),', distance from origin is ',num2str(norm(x0vec(1:2)))])
    tic
    d = x0vec(2)/x0vec(1);
    cvx_begin quiet
    variables X(n,N+1) U(m,N) Y(p,N)
    X(:,2:N+1) == Abig*X(:,1:N) + Bbig*U;
    Y(:,1:N) == Cbig*X(:,1:N) + Dbig*U;
    X(:,1) == x0vec;
    max(Y(:,1:Nc)') <= Ymax';
    max((U(1,:).^2 + U(2,:).^2)') <= Umax';
    max(U(3,:)') <= Tmax;
    min(U(3,:)') >= -Tmax;
    atan(d) - X(3,end) <= .12;
    atan(d) - X(3,end) >= -.12;
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
    Cbig = make_C(params, x0vec, 0);
    Ymax = make_Ymax(params, x0vec, 0);
    if counter>200
        disp(['More than 200 iterations. Stopping simulation.'])
        break
    end
    if params.slack_variable
        dist = sqrt((x0vec(1)-rp*cos(x0vec(3)))^2+(x0vec(2)-rp*sin(x0vec(3)))^2);
        s_weight = 10^(params.slack_slope*dist + params.slack_intercept);
        R(4,4) = s_weight; R(5,5) = s_weight;
        R(6,6) = s_weight; R(7,7) = s_weight;
        slack = [slack s_weight];
    end
    
end
% Transforming data back to inertial frame
xtotnonI = xtot(1:6,:); % Save data in non-inertial frame
xtot(1,:) = xtotnonI(1,:).*cos(xtot(8,:))-xtotnonI(2,:).*sin(xtot(8,:));
xtot(2,:) = xtotnonI(1,:).*sin(xtot(8,:))+xtotnonI(2,:).*cos(xtot(8,:));
xtot(4,:) = xtotnonI(4,:).*cos(xtot(8,:))-xtotnonI(5,:).*sin(xtot(8,:));
xtot(5,:) = xtotnonI(4,:).*sin(xtot(8,:))+xtotnonI(5,:).*cos(xtot(8,:));
xtot(3,:) = xtotnonI(3,:) + xtot(8,:);
xtot(6,:) = xtotnonI(6,:) + omega;
xtot(9,:) = abs(xtotnonI(1,:)-(rp+rs))+abs(xtotnonI(2,:));

disp('Simulation Complete')
