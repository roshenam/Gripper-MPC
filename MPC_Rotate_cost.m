function [xtot, utot,cost] = MPC_Rotate_cost(init,params,phi,omega)
% 2D accounting for platforms that rotate with a constant speed.
% Implementing cone constraints as soft through cost function.

% 6/23 Not implemented
% 6/25 Removing angle control to deal with separately later. 
% 6/27 Not functioning. Sticking with slack variable implementation for
% now.

% Extracting initial conditions and parameters from inputs
x0 = init(1); y0 = init(2); vx0 = init(3); vy0 = init(4);
rp = params.rp; rs = params.rs; 
Ts = params.Ts; 
N = params.N;
lambda = 10^2;
x0vec = [x0; y0; vx0; vy0; rp*cos(phi); rp*sin(phi); x0-rp*cos(phi);...
    y0-rp*sin(phi); phi+omega*Ts; phi];

% Creating dynamic system 
A = zeros(4,4); A(1,3) = 1; A(2,4) = 1;
B = zeros(4,2); B(3,1) = 1; B(4,2) = 1;
C = zeros(4,4);  D = zeros(4,2);
sysD = ss(A,B,C,D);
sysD = c2d(sysD,Ts);

% Adding additional states to the system to account for rotation of the
% platform. States will be rx and ry (the x and y pos of the selected
% point), and sigx and sigy (the x and y distance between the selected
% point and the center of mass of the spacecraft). Also adding states to
% predict the future state of the LOS cone constraints. State vector is now
% [x y x' y' rx ry sigx sigy z1 z2]

Abig = zeros(10,10); Abig(1:4,1:4) = sysD.a; 
Abig(5,5) = cos(omega*Ts); Abig(5,6) = -sin(omega*Ts);
Abig(6,5) = sin(omega*Ts); Abig(6,6) = cos(omega*Ts);
Abig(7,1) = 1; Abig(7,5) = -1;
Abig(8,2) = 1; Abig(8,6) = -1;
Abig(9,9) = 2; Abig(9,10) = -1;
Abig(10,9) = 1; 
Bbig = zeros(10,2); Bbig(1:4,1:2) = sysD.b;


% Y is a vector of constrained outputs. C and D are matrices defining Y's
% dependence on the states and control inputs
Dbig = zeros(2,2);
Cbig = make_C(params, x0, y0, phi,0);

% Upper bounds on thrust inputs and slack variables
Umax = params.Umax;
Ymax = [1 1]';

% Creating weighting matrices for cost function
Q = zeros(10,10); Q(1:4,1:4) = 10^4.*eye(4); R = 10^2.*eye(2); 
% Solving infinite horizon unconstrained LQR for stability enforcement
[~,P,~] = lqrd(A,B,Q(1:4,1:4),R,Ts);
Pbig = zeros(10,10);
Pbig(1:4,1:4) = P;

n = 10; m = 2; 
utot = [];
ytot = [];
xtot = x0vec;
cost = [];
counter = 0;
while norm(x0vec(1:2))>=(rp+rs)
%for i=1:20
    counter = counter + 1;
    disp(['Running optimization ',num2str(counter),', distance from origin is ',num2str(norm(x0vec(1:2)))])
    cvx_begin quiet
    variables X(n,N+1) U(m,N)
    X(:,2:N+1) == Abig*X(:,1:N) + Bbig*U;
    %Y(:,1:N) == Cbig*X(:,1:N) + Dbig*U;
    X(:,1) == x0vec;
    %max(Y(:,1:Nc)') <= Ymax';
    max((U(1,:).^2 + U(2,:).^2)') <= Umax';
    minimize (norm(Q*X(:,1:N),'fro') + norm(R*U(:,1:N),'fro') +...
        X(:,N+1)'*Pbig*X(:,N+1) + sum([lambda lambda]*(Cbig*X(:,1:N) - ones(2,N))));
    cvx_end
    u = U(:,1);
    if isnan(u(1))
        disp(['Problem became infeasible at iteration ',num2str(counter)])
        break
    end
    x0vec = Abig*x0vec+Bbig*u;
    xtot = [xtot x0vec] ;
    utot = [utot u];
    Cbig = make_C(params, x0vec(1), x0vec(2), x0vec(end),0);
    currcost = cvx_optval;
    cost = [cost; currcost];
   
end
disp('Simulation Complete')
