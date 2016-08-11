function [xtot, utot, cost, time, ytot] = MPC_Rotate_slack(init,params,phi,omega)
% 2D accounting for platforms that rotate with a constant speed.
% Implementing cone constraints as soft through use of slack variables.

% 6/23 Not implemented
% 6/25 Removing angle control to deal with separately later. 
% 6/26 Removing velocity bounding to just check if it can contact the right
% rotating point.

% Extracting initial conditions and parameters from inputs
x0 = init(1); y0 = init(2); theta0 = init(3); 
vx0 = init(4); vy0 = init(5); thetadot0 = init(6);
rp = params.rp; rs = params.rs; 
Ts = params.Ts; 
N = params.N; Nc = params.Nc; 
nu0 = theta0-phi; nu0dot = thetadot0-omega;
n = sqrt(params.mu/(params.Ro)^3);
x0vec = [x0; y0; nu0; vx0; vy0; nu0dot; phi+omega*Ts; phi];

% Creating dynamic system
A = zeros(6,6); A(1,4) = 1; A(2,5) = 1; A(3,6) = 1;
%A(4,1) = 3*n^2; A(4,5) = 2*n; A(5,4) = -2*n; %CWH Dynamics
B = zeros(6,7); B(4,1) = 1; B(5,2) = 1; B(6,3) = 1;
C = zeros(6,6);  D = zeros(6,7);
sysD = ss(A,B,C,D);
sysD = c2d(sysD,Ts);

% Adding additional states to the system to account for rotation of the
% platform. States will be rx and ry (the x and y pos of the selected
% point), and sigx and sigy (the x and y distance between the selected
% point and the center of mass of the spacecraft). Also adding states to
% predict the future state of the LOS cone constraints. State vector is now
% [x y x' y' theta theta' z1 z2 t]
Abig = zeros(8,8); Abig(1:6,1:6) = sysD.a; 
Abig(7,7) = 2; Abig(7,8) = -1; Abig(8,7) = 1;
Bbig = zeros(8,7); Bbig(1:6,1:3) = sysD.b(1:6,1:3);


% Y is a vector of constrained outputs. C and D are matrices defining Y's
% dependence on the states and control inputs

%Dbig = zeros(4,7); Dbig(1,4) = 1; Dbig(2,5) = 1; % Slack variables to make constraints soft
%Dbig(3,6) = 1; Dbig(4,7) = 1;
Dbig = zeros(4,7); Dbig(1,4) = 1; Dbig(2,5) = 1; 
%Dbig(3,6) = 1; Dbig(4,7) = 1;
Cbig = make_C(params, x0vec, 1);


% Upper bounds on thrust inputs and slack variables
Umax = params.Umax;
Tmax = params.Tmax; % max torque applied by flywheel
Ymax = make_Ymax(params, x0vec,1);

% Creating weighting matrices for cost function
Q = params.Qval.*eye(8); Q(7,7) = 0; Q(8,8) = 0; R = params.Rval.*eye(7); 
R(4,4) = params.slackweight; R(5,5) = params.slackweight;
R(6,6) = params.slackweight; R(7,7) = params.slackweight;
% Solving infinite horizon unconstrained LQR for stability enforcement
[~,P,~] = lqrd(A,B(1:6,1:3),Q(1:6,1:6),R(1:3,1:3),Ts);
Pbig = zeros(8,8);
Pbig(1:6,1:6) = P;

n = 8; m = 7; p=4; 
time = [];
utot = [];
ytot = [Ymax];
xtot = x0vec;
cost = [];
counter = 0;
%Xr = zeros(2,N+1);
%Xr(2,:) = ones.*omega;
%for k=1:N+1
%Xr(1,k) = phi + k*Ts*omega;
%end
%Mtrack = [10000 0; 100 0];
lambda1 = [100 100 100 10^5 10^5 10^5 10^5];
lambda2 = 100.*ones(1,6);
while norm(x0vec(1:2))>=(rp+rs)
%for j=1:50
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
    max(U(3,:)') <= Tmax;
    min(U(3,:)') >= -Tmax;
    %if (phi-theta0)>0
    %max(X(7,:)') <= x0vec(6)+N*Ts*omega;
    %else
    %    min(X(7,:)') >= x0vec(6)-N*Ts*omega;
    %end
    minimize (norm(Q*X(:,1:N),'fro') + norm(R*U(:,1:N),'fro') +...
        X(:,N+1)'*Pbig*X(:,N+1)) %+ norm(Mtrack*(X(7:8,:)-Xr),'fro'));
        %minimize (sum(lambda1*abs(U(:,1:N))) + sum(lambda2*abs(X(1:6,1:N))))
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
    Cbig = make_C(params, x0vec, 1);
    Ymax = make_Ymax(params, x0vec, 1);
    ytot = [ytot Ymax];
    currcost = cvx_optval;
    cost = [cost; currcost];
     if counter>100
        disp(['More than 200 iterations. Stopping simulation.'])
        break
    end
    %Xr = Xr + Ts*omega;
   
end
xtot(3,:) = xtot(3,:) + xtot(8,:);
xtot(6,:) = xtot(6,:) + omega;
xtot(9,:) = abs(xtot(1,:)-(rp+rs).*cos(xtot(8,:)))+abs(xtot(2,:)-...
    (rp+rs).*sin(xtot(8,:)));
xtot(10,:) = -(xtot(4,:).*cos(xtot(8,:))+xtot(5,:).*sin(xtot(8,:)));

disp('Simulation Complete')
