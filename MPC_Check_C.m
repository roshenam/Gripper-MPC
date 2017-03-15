function [xtot, utot, cost, time, phidata] = MPC_Check_C(init,params)
% 2D accounting for platforms that rotate with a constant speed.
% Implementing cone constraints as soft through use of slack variables.

% 6/23 Not implemented
% 6/25 Removing angle control to deal with separately later.
% 6/26 Removing velocity bounding to just check if it can contact the right
% rotating point.
% 8/24 Making phi not a state, rewriting cone constraints so they're not
% approximations

% Extracting initial conditions and parameters from inputs
rtol = params.rtol;
phi = params.phi; omega = params.omega;
gamma = params.gamma;
x0 = init(1); y0 = init(2); theta0 = init(3);
vx0 = init(4); vy0 = init(5); thetadot0 = init(6);
rp = params.rp; rs = params.rs;
Ts = params.Ts;
N = params.N; Nc = params.Nc;
nu0 = theta0-phi; nu0dot = thetadot0-omega;
%n = sqrt(params.mu/(params.Ro)^3); % Parameter for CWH dynamics
%x0vec = [x0; y0; nu0; vx0; vy0; nu0dot; phi+omega*Ts; phi];
x0vec = [x0; y0; nu0; vx0; vy0; nu0dot];

A = params.A; B = params.B;

% Upper bounds on thrust inputs and slack variables
Umax = params.Umax;
Tmax = params.Tmax; % max torque a=pplied by flywheel
%Ymax = make_Ymax(params, x0vec,1);

% Creating weighting matrices for cost function
Q = params.Q;
R = params.R;
P = params.P;
cA, cb] = Calc_InvSet(params, Abig, Bbig(:,1:3));

n = 6; m = 5; %p=4;
time = [];
utot = [];
%ytot = [Ymax];
xtot = x0vec;
cost = [];
counter = 0;
phidata = phi;
%Xr = zeros(2,N+1);
%Xr(2,:) = ones.*omega;
%for k=1:N+1
%Xr(1,k) = phi + k*Ts*omega;
%end
%Mtrack = [10000 0; 100 0];
%lambda1 = [100 100 100 10^5 10^5 10^5 10^5];
%lambda2 = 100.*ones(1,6);
multiplier = 0:1:N-1;
phis = phi + omega.*Ts.*multiplier;
x1dot_min = -0.5;       % m/s (true limit is ~0.3)
x1dot_max = 0.5; 
x2dot_min = -0.5;
x2dot_max = 0.5;
% x1_min = 0.1;         % meters (true limits)
% x1_max = 3.5;   
% x2_min = 0.1;
% x2_max = 2.5;
x1_min = -10;           % meters (relaxed limits)
x1_max = 10;   
x2_min = -10;
x2_max = 10;
thetadot_min = -pi;     % rad/s (made up)
thetadot_max = pi; 
theta_min = -50*pi;     % rad (in principle, should have no limits)
theta_max = 50*pi;
xmax = [x1_max x2_max theta_max x1dot_max x2dot_max thetadot_max]';
xmin = [x1_min x2_min theta_min x1dot_min x2dot_min thetadot_min]';

%while norm(x0vec(1:2))>=(rp+rs)
    for j=1
    
    counter = counter + 1;
    disp(['Running optimization ',num2str(counter),', distance from origin is ',num2str(norm(x0vec(1:2)))])
    tic
    cvx_begin quiet
    variables X(n,N+1) U(m,N)
    X(:,2:N+1) == Abig*X(:,1:N) + Bbig*U;
    %Y(:,1:N) == Cbig*X(:,1:N) + Dbig*U;
    X(:,1) == x0vec;
    %max(Y(:,1:Nc)') <= Ymax';
    for j = 1:N
        L1 = sin(phis(j)+gamma)./((rp-rtol)*sin(gamma));
        L2 = -cos(phis(j)+gamma)./((rp-rtol)*sin(gamma));
        L3 = sin(phis(j)-gamma)./((rp-rtol)*sin(gamma));
        L4 = -cos(phis(j)-gamma)./((rp-rtol)*sin(gamma));
        X(1,j)*L1 + X(2,j)*L2 + U(4,j) <= 1;
        X(1,j)*L3 + X(2,j)*L4 + U(5,j) >= 1; 
    end
    max(U(1,:)) <= Umax;
    min(U(1,:)) >= -Umax;
    max(U(2,:)) <= Umax;
    min(U(2,:)) >= -Umax;
    max(U(3,:)') <= Tmax;
    min(U(3,:)') >= -Tmax;
    max(X') <= xmax';
    min(X') >= xmin';
    max((cA*X(:,N+1))') <= cb';
    minimize (norm(Q*X(:,1:N),'fro') + norm(R*U(:,1:N),'fro') +...
        X(:,N+1)'*Pbig*X(:,N+1)) %+ norm(Mtrack*(X(7:8,:)-Xr),'fro'));
    %minimize (sum(lambda1*abs(U(:,1:N))) + sum(lambda2*abs(X(1:6,1:N))))
    cvx_end
    if strcmp(cvx_status,'Infeasible')
        %time = 0:1:counter;
        %phidata = phi + omega.*Ts.*time;
        %xtot(3,:) = xtot(3,:) + phidata;
        %xtot(6,:) = xtot(6,:) + omega;
        %xtot(7,:) = phidata;
        disp(['Problem became infeasible'])
        return
    end
    timecurr = toc;
    u = U(:,1);
    if isnan(u(1))
        disp(['Problem became infeasible at iteration ',num2str(counter)])
        break
    end
    phidata = [phidata phis(2)];
    x0vec = Abig*x0vec+Bbig*u;
    phis = phis + omega.*Ts;
    time = [time timecurr];
    xtot = [xtot x0vec] ;
    utot = [utot u];
    %Cbig = make_C(params, x0vec, 1);
    %Ymax = make_Ymax(params, x0vec, 1);
    %ytot = [ytot Ymax];
    currcost = cvx_optval;
    cost = [cost; currcost];
    if counter>100
        disp(['More than 200 iterations. Stopping simulation.'])
        break
    end
    %Xr = Xr + Ts*omega;
    
end
time = 0:1:counter;
phidata = phi + omega.*Ts.*time;
xtot(3,:) = xtot(3,:) + phidata;
xtot(6,:) = xtot(6,:) + omega;
xtot(7,:) = phidata;
%xtot(9,:) = abs(xtot(1,:)-(rp+rs).*cos(xtot(8,:)))+abs(xtot(2,:)-...
%(rp+rs).*sin(xtot(8,:)));
%xtot(10,:) = -(xtot(4,:).*cos(xtot(8,:))+xtot(5,:).*sin(xtot(8,:)));

disp('Simulation Complete')
