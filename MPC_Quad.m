function [xtot, utot,ytot] = MPC_Quad(init,params,phi)
% This version uses a quadratic cone instead of a linear cone. Turns out
% this is not convex. Will explore in the future. For now have left in the
% non-convex quadratic constraints and have added quadratic constraints for
% the thrust max
% Extracting initial conditions and parameters from inputs
x0 = init(1); y0 = init(2); theta0 = init(3); vx0 = init(4); vy0 = init(5);
omega0 = init(6); 
rp = params.rp; rtol = params.rtol; rs = params.rs; phioff = params.gamma;
Ts = params.Ts; eta = params.eta; etaR = params.etaR; 
beta1 = params.beta1; beta2 = params.beta2; beta3 = params.beta3; 
beta4 = params.beta4; beta5 = params.beta5; 
N = params.N; Nc = params.Nc; 

x0vec = [x0; y0; theta0-phi; vx0; vy0; omega0];

% Computing quadratic cone constants
aA = phi; aB = aA + pi;
l = rtol-rp; 
xd1 = rp*cos(phi+phioff); yd1 = rp*sin(phi+phioff);
xd2 = rp*cos(phi-phioff); yd2 = rp*sin(phi-phioff);
kA = (yd1*cos(aA)-xd1*sin(aA)-l)/(xd1*cos(aA)+yd1*sin(aA))^2;
kB = (yd2*cos(aB)-xd2*sin(aB)-l)/(xd2*cos(aB)+yd2*sin(aB))^2;

% Creating dynamic system 
A = zeros(6,6); A(1,4) = 1; A(2,5) = 1; A(3,6) = 1; 
B = zeros(6,9); B(4,1) = 1; B(5,2) = 1; B(6,3) = 1;

% Y is a vector of linearly constrained outputs, which here are only due to
% the velocity constraints. The cone constraints are now quadratic and are
% not linear. C and D are matrices defining Y's dependence on the states 
% and control inputs
C = zeros(6,6);  D = zeros(6,9);

C(1,4) = eta*cos(theta0); C(1,5) = eta*sin(theta0); % Vn bounds
C(2,4) = -eta*cos(theta0); C(2,5) = -eta*sin(theta0);
C(3,4) = -eta*sin(theta0); C(3,5) = eta*cos(theta0); % Vt bounds
C(4,4) = eta*sin(theta0); C(4,5) = -eta*cos(theta0);
C(5,6) = eta; C(6,6) = -eta; % Theta' bounds

D(1,4) = eta; D(2,5) = eta; D(3,6) = eta; D(4,7) = eta; 
D(5,8) = eta; D(6,9) = eta; 

% Upper bounds on thrust inputs and slack variables
Umax = [0.2 0.2 0.2 10^10 10^10 10^10 10^10 10^10 10^10]';
Umin = -Umax;
maxterm = abs(x0-(rp+rs)*cos(phi))+abs(y0-(rp+rs)*sin(phi));
Ymax = [maxterm+beta1 maxterm+beta2 maxterm+beta3 maxterm+beta3...
    etaR*maxterm+beta4 etaR*maxterm+beta5]';
sysD = ss(A,B,C,D);
sysD = c2d(sysD,Ts);
% Creating weighting matrices for cost function
Q = 10^3.*eye(6); Q(3,3) = 100; R = 10^6.*eye(9); 
R(1,1) = 100; R(2,2) = 100; R(3,3) = 100;
% Solving infinite horizon unconstrained LQR for stability enforcement
[K,P,~] = lqrd(A,B(:,1:3),Q,R(1:3,1:3),Ts);

n = 6; m = 9; p = 6; 
utot = [];
ytot = [];
xtot = x0vec;
counter = 0;
% changes that don't matter
while norm(x0vec(1:2))>=(rp+rs)
    counter = counter + 1;
    disp(['Running optimization ',num2str(counter),', distance from origin is ',num2str(norm(x0vec(1:2)))])
    cvx_begin quiet
    variables X(n,N+1) U(m,N) Y(p,N)
    X(:,2:N+1) == sysD.a*X(:,1:N) + sysD.b*U;
    Y(:,1:N) == sysD.c*X(:,1:N) + sysD.d*U;
    X(:,1) == x0vec;
    max(Y(:,1:Nc)') <= Ymax';
    max(U(1,:).^2+U(2,:).^2) <= Umax(1);
    min(-(U(1,:).^2+U(2,:).^2)) >= -Umax(1);
    max(U(3,:)) <= Umax(3);
    min(-U(3,:)) >= -Umax(3);
    %min(U') >= Umin';
    %max(U') <= Umax';
    min(kA.*(X(1,1:Nc).*cos(aA)+X(2,1:Nc).*sin(aA)).^2+l+X(1,1:Nc).*sin(aA)-X(2,1:Nc).*cos(aA)) >= 0;
    min(kB.*(X(1,1:Nc).*cos(aB)+X(2,1:Nc).*sin(aB)).^2+l+X(1,1:Nc).*sin(aB)-X(2,1:Nc).*cos(aB)) >= 0;
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
    C(1,4) = eta*cos(x0vec(3)+phi); C(1,5) = eta*sin(x0vec(3)+phi);
    C(2,4) = -eta*cos(x0vec(3)+phi); C(2,5) = -eta*sin(x0vec(3)+phi);
    C(3,4) = -eta*sin(x0vec(3)+phi); C(3,5) = eta*cos(x0vec(3)+phi);
    C(4,4) = eta*sin(x0vec(3)+phi); C(4,5) = -eta*cos(x0vec(3)+phi);
    maxtermcurr = abs(x0vec(1)-(rp+rs)*cos(phi))+abs(x0vec(2)-(rp+rs)*sin(phi));
    Ymax = [maxtermcurr+beta1 maxtermcurr+beta2 maxtermcurr+beta3 ...
        maxtermcurr+beta3 etaR*maxtermcurr+beta4 etaR*maxtermcurr+beta5]';
   
end
disp('Simulation Complete')
