function [xout, xtotnonI, utot, cost, time, ytot] = MPC_nonI_invset(init,params,phi,omega)
% 2D accounting for platforms that rotate with a constant speed.
% Implementing cone constraints as soft through use of slack variables.
% This version is modeled in the non-inertial reference frame to make
% constraints constant

% Extracting initial conditions and parameters from inputs
x0 = init(1); y0 = init(2); theta0 = init(3); 
vx0 = init(4); vy0 = init(5); thetadot0 = init(6);
rp = params.rp; rs = params.rs; 
vD = params.vD;
Ts = params.Ts; 
N = params.N; Nc = params.Nc; 
gamma = params.gamma;
phi = params.phi;
omega = params.omega;
Rmat = [cos(phi) sin(phi); -sin(phi) cos(phi)];
r0 = Rmat*[x0;y0];
v0 = Rmat*[vx0;vy0];
nu0 = theta0-phi; nu0dot = thetadot0-omega;
% 7th state is time, 8th state is dt/dt, 9th state is phi(j+1), 10th state
% is phi(j)
x0vec = [r0(1); r0(2); nu0; v0(1)+vD; v0(2); nu0dot; 0; 1; phi+omega*Ts; phi];

% Creating dynamic system
A = zeros(8,8); A(1,4) = 1; A(2,5) = 1; A(3,6) = 1;
A(4,1) = omega^2; A(4,5) = 2*omega; 
A(5,2) = omega^2; A(5,4) = -2*omega;
A(7,8) = 1; % 7th state is time, 8th state is dt/dt
B = zeros(8,3); B(4,1) = 1; B(5,2) = 1; B(6,3) = 1;
C = zeros(8,8);  D = zeros(8,3);
sysD = ss(A,B,C,D);
sysD = c2d(sysD,Ts);

% Adding additional states to the system to account for rotation of the
% platform.
Abig = zeros(10,10); Abig(1:8,1:8) = sysD.a; 
Abig(9,9) = 2; Abig(9,10) = -1; Abig(10,9) = 1;
Bbig = zeros(10,3); Bbig(1:6,1:3) = sysD.b(1:6,1:3);

% Y is a vector of constrained outputs. C and D are matrices defining Y's
% dependence on the states and control inputs

Dbig = zeros(2,3);
Cbig = make_C(params, x0vec, 2);

% Upper bounds on thrust inputs and slack variables
Umax = params.Umax;
Tmax = params.Tmax; % max torque applied by flywheel
Ymax = make_Ymax(params, x0vec, 2);

% Creating weighting matrices for cost function
Q = params.Qval.*eye(10); Q(7,7) = 0; Q(8,8) = 0; Q(9,9) = 0; Q(10,10) = 0; 
R = params.Rval.*eye(3);

% Solving infinite horizon unconstrained LQR for stability enforcement
[~,P,~] = lqrd(A(1:6,1:6),B(1:6,1:3),Q(1:6,1:6),R(1:3,1:3),Ts);
Pbig = zeros(10,10);
Pbig(1:6,1:6) = P;

n = 10; m = 3; p=2; 
time = [];
utot = [];
ytot = Ymax;
xtot = x0vec;
cost = [];
counter = 0;
while norm([x0vec(1)-vD*counter*Ts,x0vec(2)])>=(rp+rs)
%for j=1:50
    counter = counter + 1;
    disp(['Running optimization ',num2str(counter),', distance from origin is ',num2str(norm([x0vec(1)-vD*counter*Ts,x0vec(2)]))])
    tic
    cvx_begin quiet
    variables X(n,N+1) U(m,N) Y(p,N)
    X(:,2:N+1) == Abig*X(:,1:N) + Bbig*U;
    Y(:,1:N) == Cbig*X(:,1:N) + Dbig*U;
    X(:,1) == x0vec;
    max(Y') <= Ymax';
    max((U(1,:).^2 + U(2,:).^2)') <= Umax';
    max(U(3,:)') <= Tmax;
    min(U(3,:)') >= -Tmax;
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
    if counter>100
        disp(['More than 100 iterations. Stopping simulation.'])
        break
    end
    
end
% Transforming data back to inertial frame
xtotnonI = xtot(1:6,:); % Save data in non-inertial frame
xreal = xtot(1,:) - vD.*xtot(7,:);
vxreal = xtot(4,:) - vD;
[dim,num] = size(xtot);
xout = zeros(7,num);
xout(1,:) = xreal.*cos(xtot(10,:))-xtot(2,:).*sin(xtot(10,:));
xout(2,:) = xreal.*sin(xtot(10,:))+xtot(2,:).*cos(xtot(10,:));
xout(3,:) = xtot(3,:) + xtot(10,:);
xout(4,:) = vxreal.*cos(xtot(10,:))-xtot(5,:).*sin(xtot(10,:));
xout(5,:) = vxreal.*sin(xtot(10,:))+xtot(5,:).*cos(xtot(10,:));
xout(6,:) = xtot(6,:) + omega;
xout(7,:) = xtot(end,:);
xout(8,:) = xtot(7,:);
%xtot(9,:) = abs(xtotnonI(1,:)-(rp+rs))+abs(xtotnonI(2,:));

disp('Simulation Complete')
