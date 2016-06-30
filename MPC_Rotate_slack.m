function [xtot, utot, cost, time] = MPC_Rotate_slack(init,params,phi,omega)
% 2D accounting for platforms that rotate with a constant speed.
% Implementing cone constraints as soft through use of slack variables.

% 6/23 Not implemented
% 6/25 Removing angle control to deal with separately later. 
% 6/26 Removing velocity bounding to just check if it can contact the right
% rotating point.

% Extracting initial conditions and parameters from inputs
x0 = init(1); y0 = init(2); vx0 = init(3); vy0 = init(4);
rp = params.rp; rs = params.rs; 
Ts = params.Ts; 
N = params.N; Nc = params.Nc; 

x0vec = [x0; y0; vx0; vy0; phi+omega*Ts; phi];

% Creating dynamic system 
A = zeros(4,4); A(1,3) = 1; A(2,4) = 1;
B = zeros(4,4); B(3,1) = 1; B(4,2) = 1;
C = zeros(4,4);  D = zeros(4,4);
sysD = ss(A,B,C,D);
sysD = c2d(sysD,Ts);

% Adding additional states to the system to account for rotation of the
% platform. States will be rx and ry (the x and y pos of the selected
% point), and sigx and sigy (the x and y distance between the selected
% point and the center of mass of the spacecraft). Also adding states to
% predict the future state of the LOS cone constraints. State vector is now
% [x y x' y' rx ry sigx sigy z1 z2]

Abig = zeros(6,6); Abig(1:4,1:4) = sysD.a; 
Abig(5,5) = 2; Abig(5,6) = -1;
Abig(6,5) = 1; 
Bbig = zeros(6,4); Bbig(1:4,1:4) = sysD.b;


% Y is a vector of constrained outputs. C and D are matrices defining Y's
% dependence on the states and control inputs
Dbig = zeros(2,4); Dbig(1,3) = 1; Dbig(2,4) = 1; % Slack variables to make constraints soft
Cbig = make_C(params, x0, y0, phi,1);

% Upper bounds on thrust inputs and slack variables
Umax = params.Umax;
Ymax = make_Ymax(params, x0, y0, phi);

% Creating weighting matrices for cost function
Q = zeros(6,6); Q(1:4,1:4) = params.Qval.*eye(4); R = params.Rval.*eye(4); 
R(3,3) = params.slackweight; R(4,4) = params.slackweight;
% Solving infinite horizon unconstrained LQR for stability enforcement
[~,P,~] = lqrd(A,B(:,1:2),Q(1:4,1:4),R(1:2,1:2),Ts);
Pbig = zeros(6,6);
Pbig(1:4,1:4) = P;

n = 6; m = 4; p=2; 
time = [];
utot = [];
ytot = [];
xtot = x0vec;
cost = [];
counter = 0;
while norm(x0vec(1:2))>=(rp+rs)
%for j=1:100
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
    minimize (norm(Q*X(:,1:N),'fro') + norm(R*U(:,1:N),'fro') + X(:,N+1)'*Pbig*X(:,N+1));
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
    Cbig = make_C(params, x0vec(1), x0vec(2), x0vec(end),1);
    Ymax = make_Ymax(params, x0vec(1), x0vec(2), x0vec(end));
    currcost = cvx_optval;
    cost = [cost; currcost];
   
end
disp('Simulation Complete')
