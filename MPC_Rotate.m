function [xtot, utot,ytot] = MPC_Rotate(init,params,phi,omega)
% 2D accounting for platforms that rotate with a constant speed. 

% 6/23 Not implemented
% 6/25 Removing angle control to deal with separately later. 
% 6/26 Removing velocity bounding to just check if it can contact the right
% rotating point.

% Extracting initial conditions and parameters from inputs
x0 = init(1); y0 = init(2); vx0 = init(3); vy0 = init(4);
rp = params.rp; rs = params.rs; 
Ts = params.Ts; 
N = params.N; Nc = params.Nc; 

x0vec = [x0; y0; vx0; vy0; rp*cos(phi); rp*sin(phi); x0-rp*cos(phi);...
    y0-rp*sin(phi); phi+omega*Ts; phi];

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

Abig = zeros(10,10); Abig(1:4,1:4) = sysD.a; 
Abig(5,5) = cos(omega*Ts); Abig(5,6) = -sin(omega*Ts);
Abig(6,5) = sin(omega*Ts); Abig(6,6) = cos(omega*Ts);
Abig(7,1) = 1; Abig(7,5) = -1;
Abig(8,2) = 1; Abig(8,6) = -1;
Abig(9,9) = 2; Abig(9,10) = -1;
Abig(10,9) = 1; 
Bbig = zeros(10,4); Bbig(1:4,1:4) = sysD.b;


% Y is a vector of constrained outputs. C and D are matrices defining Y's
% dependence on the states and control inputs
Dbig = zeros(2,4); Dbig(1,3) = 1; Dbig(2,4) = 1; % Slack variables to make constraints soft
Cbig = make_C(params, x0, y0, phi);

% Upper bounds on thrust inputs and slack variables
Umax = params.Umax;
Ymax = [1 1]';

% Creating weighting matrices for cost function
Q = zeros(10,10); Q(1:4,1:4) = 10^3.*eye(4); R = 10^2.*eye(4); 
R(3,3) = 10^4; R(4,4) = 10^4;
% Solving infinite horizon unconstrained LQR for stability enforcement
[~,P,~] = lqrd(A,B(:,1:2),Q(1:4,1:4),R(1:2,1:2),Ts);
Pbig = zeros(10,10);
Pbig(1:4,1:4) = P;

n = 10; m = 4; p=2; 
utot = [];
ytot = [];
xtot = x0vec;
counter = 0;
while norm(x0vec(1:2))>=(rp+rs)
    counter = counter + 1;
    disp(['Running optimization ',num2str(counter),', distance from origin is ',num2str(norm(x0vec(1:2)))])
    cvx_begin quiet
    variables X(n,N+1) U(m,N) Y(p,N)
    X(:,2:N+1) == Abig*X(:,1:N) + Bbig*U;
    Y(:,1:N) == Cbig*X(:,1:N) + Dbig*U;
    X(:,1) == x0vec;
    max(Y(:,1:Nc)') <= Ymax';
    max((U(1,:).^2 + U(2,:).^2)') <= Umax';
    minimize (norm(Q*X(:,1:N),'fro') + norm(R*U(:,1:N),'fro') + X(:,N+1)'*Pbig*X(:,N+1));
    cvx_end
    u = U(:,1);
    if isnan(u(1))
        disp(['Problem became infeasible at iteration ',num2str(counter)])
        break
    end
    x0vec = Abig*x0vec+Bbig*u;
    xtot = [xtot x0vec] ;
    utot = [utot u];
    Cbig = make_C(params, x0vec(1), x0vec(2), x0vec(end));
    
   
end
disp('Simulation Complete')
