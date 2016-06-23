function [xtot, utot,ytot] = MPC_Rotate(init,params,phi,omega,pred)
% 2D accounting for platforms that rotate with a constant speed. 

% 6/23 Not implemented

% If pred is zero, there will be no prediction of the LOS cone constraints
% in the model. If pre is one, there will be prediction of the LOS cone
% constraints

% Extracting initial conditions and parameters from inputs
x0 = init(1); y0 = init(2); theta0 = init(3); vx0 = init(4); vy0 = init(5);
omega0 = init(6); 
rp = params.rp; rtol = params.rtol; rs = params.rs; gamma = params.gamma;
Ts = params.Ts; eta = params.eta; etaR = params.etaR; 
beta1 = params.beta1; beta2 = params.beta2; beta3 = params.beta3; 
beta4 = params.beta4; beta5 = params.beta5; 
N = params.N; Nc = params.Nc; 

x0vec = [x0; y0; theta0-phi; vx0; vy0; omega0];

% Creating dynamic system 
A = zeros(6,6); A(1,4) = 1; A(2,5) = 1; A(3,6) = 1; 
B = zeros(6,11); B(4,1) = 1; B(5,2) = 1; B(6,3) = 1;

% Adding additional states to the system to account for rotation of the
% platform. States will be rx and ry (the x and y pos of the selected
% point), and sigx and sigy (the x and y distance between the 
% Y is a vector of constrained outputs. C and D are matrices defining Y's
% dependence on the states and control inputs
C = zeros(8,6);  D = zeros(8,11);


C(1,1) = sin(phi+gamma)/((rp-rtol)*sin(gamma)); % LOS cone constraints
C(1,2) = -cos(phi+gamma)/((rp-rtol)*sin(gamma));
C(2,1) = -sin(phi-gamma)/((rp-rtol)*sin(gamma));
C(2,2) = cos(phi-gamma)/((rp-rtol)*sin(gamma));
C(3,4) = eta*cos(theta0); C(3,5) = eta*sin(theta0); % Vn bounds
C(4,4) = -eta*cos(theta0); C(4,5) = -eta*sin(theta0);
C(5,4) = -eta*sin(theta0); C(5,5) = eta*cos(theta0); % Vt bounds
C(6,4) = eta*sin(theta0); C(6,5) = -eta*cos(theta0);
C(7,6) = eta; C(8,6) = -eta; % Theta' bounds

D(1,4) = eta; D(2,5) = eta; D(3,6) = eta; D(4,7) = eta; 
D(5,8) = eta; D(6,9) = eta; D(7,10) = eta; D(8,11) = eta;

% Upper bounds on thrust inputs and slack variables
Umax = [0.2 0.2 0.2 10^10 10^10 10^10 10^10 10^10 10^10 10^10 10^10]';
Umin = -Umax;
maxterm = abs(x0-(rp+rs)*cos(phi))+abs(y0-(rp+rs)*sin(phi));
Ymax = [1 1 maxterm+beta1 maxterm+beta2 maxterm+beta3 maxterm+beta3...
    etaR*maxterm+beta4 etaR*maxterm+beta5]';
sysD = ss(A,B,C,D);
sysD = c2d(sysD,Ts);
% Creating weighting matrices for cost function
Q = 10^3.*eye(6); Q(3,3) = 100; R = 10^6.*eye(11); 
R(1,1) = 100; R(2,2) = 100; R(3,3) = 100;
% Solving infinite horizon unconstrained LQR for stability enforcement
[K,P,~] = lqrd(A,B(:,1:3),Q,R(1:3,1:3),Ts);

n = 6; m = 11; p = 8; 
utot = [];
ytot = [];
xtot = x0vec;
counter = 0;
while norm(x0vec(1:2))>=(rp+rs)
    counter = counter + 1;
    disp(['Running optimization ',num2str(counter),', distance from origin is ',num2str(norm(x0vec(1:2)))])
    cvx_begin quiet
    variables X(n,N+1) U(m,N) Y(p,N)
    X(:,2:N+1) == sysD.a*X(:,1:N) + sysD.b*U;
    Y(:,1:N) == sysD.c*X(:,1:N) + sysD.d*U;
    X(:,1) == x0vec;
    max(Y(:,1:Nc)') <= Ymax';
    min(U') >= Umin';
    max(U') <= Umax';
    minimize (norm(Q*X(:,1:N),'fro') + norm(R*U(:,1:N),'fro') + X(:,N+1)'*P*X(:,N+1));
    cvx_end
    u = U(:,1);
    if isnan(u(1))
        disp(['Problem became infeasible at iteration ',num2str(counter)])
        break
    end
    %y = Y(:,1);
    x0vec = sysD.a*x0vec+sysD.b*u;
    xtot = [xtot x0vec] ;
    utot = [utot u];
    %ytot = Y;
    % Updating constraints for velocity bounding
    C(3,4) = eta*cos(x0vec(3)+phi); C(3,5) = eta*sin(x0vec(3)+phi);
    C(4,4) = -eta*cos(x0vec(3)+phi); C(4,5) = -eta*sin(x0vec(3)+phi);
    C(5,4) = -eta*sin(x0vec(3)+phi); C(5,5) = eta*cos(x0vec(3)+phi);
    C(6,4) = eta*sin(x0vec(3)+phi); C(6,5) = -eta*cos(x0vec(3)+phi);
    maxtermcurr = abs(x0vec(1)-(rp+rs)*cos(phi))+abs(x0vec(2)-(rp+rs)*sin(phi));
    Ymax = [1 1 maxtermcurr+beta1 maxtermcurr+beta2 maxtermcurr+beta3 ...
        maxtermcurr+beta3 etaR*maxtermcurr+beta4 etaR*maxtermcurr+beta5]';
   
end
disp('Simulation Complete')
