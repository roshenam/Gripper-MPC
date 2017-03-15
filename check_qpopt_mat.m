
folder_name = 'FF_filegen/';
fname_param = 'FF_params_test.h';
fname_mpc = 'FF_test.dat';

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

% Creating dynamic system
% States are x, y, theta, xdot, ydot, thetadot
A = zeros(6,6); A(1,4) = 1; A(2,5) = 1; A(3,6) = 1;
%A(4,1) = 3*n^2; A(4,5) = 2*n; A(5,4) = -2*n; %CWH Dynamics
B = zeros(6,5); B(4,1) = 1; B(5,2) = 1; B(6,3) = 1;
C = zeros(6,6);  D = zeros(6,5);
sysD = ss(A,B,C,D);
sysD = c2d(sysD,Ts);
Ad = sysD.a; Bd = sysD.b;

% Upper bounds on thrust inputs and slack variables
Umax = params.Umax;
Tmax = params.Tmax; % max torque a=pplied by flywheel
%Ymax = make_Ymax(params, x0vec,1);

% Creating weighting matrices for cost function
Q = params.Qval.*eye(6); %Q(7,7) = 0; Q(8,8) = 0;
R = params.Rval.*eye(5);
R(4,4) = params.slackweight; R(5,5) = params.slackweight;
%R(6,6) = params.slackweight; R(7,7) = params.slackweight;
% Solving infinite horizon unconstrained LQR for stability enforcement
[~,P,~] = lqrd(A,B(1:6,1:3),Q(1:6,1:6),R(1:3,1:3),Ts);
%Pbig = zeros(8,8);
Pbig(1:6,1:6) = P;

n = 6; m = 5;
multiplier = 0:1:N-1;
phis = phi + omega.*Ts.*multiplier;

% Concatenates all constraints due to dynamics, thrust, state limits, and
% terminal conditions into one matrix and vector. Initial conditions and
% cone constraints will be tacked on in the C code. Will stack all
% variables on top of each other to create one vector of parameters. Order
% will be [x[0]; x[1]; x[2]; ... x[N]; u[0]; u[1]; ... u[N-1]]. With six
% states, three inputs, and N time steps, we'll have 6*(N+1)+3*(N)
% variables. x is ordered as [xdot ydot x y theta thetadot] and u is
% ordered as [fx fy tz]. l and u are upper and lower bounds 
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
u1_min = -.2;            % Two thrusters, normalized
u1_max = .2;             % Cut actuation limit in half due to 50% actuation duty cycle
u2_min = -.2;
u2_max = .2;
u3_min = -.5;    % [rad/s] Desired rotation rate
u3_max = .5;

% Set final state constraint to maximally invariant set
model = LTISystem('A',Ad,'B',Bd(:,1:3));
model.x.min = [x1_min x2_min theta_min x1dot_min x2dot_min thetadot_min]';
model.x.max = [x1_max x2_max theta_max x1dot_max x2dot_max thetadot_max]';
model.u.min = [u1_min u2_min u3_min];
model.u.max = [u1_max u2_max u3_max];
invset = model.invariantSet('maxIterations',40);
cA = invset.A; cb = invset.b;

rowdim = n*N; coldim = n*(N+1)+m*N;
% Constructing matrix of all dynamics constraints
Adyn = zeros(rowdim,coldim);

for t = 1:N
    Adyn(n*(t-1)+1:n*(t-1)+n, n*(t-1)+1:n*(t-1)+n) = Ad;
    Adyn(n*(t-1)+1:n*(t-1)+n, n*t+1:n*t+n) = -eye(n);
    Adyn(n*(t-1)+1:n*(t-1)+n, n*(N+1)+m*(t-1)+1:n*(N+1)+m*t) = Bd;
end

% Constructing vector of all dynamics constraints
bdyn = zeros(rowdim,1);
ldyn = bdyn; udyn = bdyn;

% Constructing vector of all state and thrust constraints
lXlims = zeros(n*(N+1),1);
uXlims = zeros(n*(N+1),1);
lUlims = zeros(m*N,1);
uUlims = zeros(m*N,1);

h_u = [u1_max u2_max u3_max 10^5 10^5 u1_min u2_min u3_min -10^5 -10^5]';
h_x = [x1_max x2_max theta_max x1dot_max  x2dot_max thetadot_max ...
      x1_min x2_min theta_min x1dot_min x2dot_min thetadot_min]';
  
for i=1:N+1
    uXlims(n*(i-1)+1:n*i,1) = h_x(1:n);
    lXlims(n*(i-1)+1:n*i,1) = h_x(n+1:2*n);
    if i<=N
        uUlims(m*(i-1)+1:m*i,1) = h_u(1:m);
        lUlims(m*(i-1)+1:m*i,1) = h_u(m+1:2*m);
    end
end
llims = [lXlims; lUlims];
ulims = [uXlims; uUlims];

% Constructing matrix of terminal constraints
nTc = length(cb);
Aterm = zeros(nTc, coldim);
Aterm(:, n*N+1:n*N+n) = cA;

% Constructing vector of terminal constraints
uterm = cb;

% Matrix for LOS cone constraints
ALOS = zeros(2*N, coldim);
for j = 1:N
        L1 = sin(phis(j)+gamma)./((rp-rtol)*sin(gamma));
        L2 = -cos(phis(j)+gamma)./((rp-rtol)*sin(gamma));
        L3 = sin(phis(j)-gamma)./((rp-rtol)*sin(gamma));
        L4 = -cos(phis(j)-gamma)./((rp-rtol)*sin(gamma));
        ALOS(2*(j-1)+1, n*(j-1)+1) = L1;
        ALOS(2*(j-1)+1, n*(j-1)+2) = L2;
        ALOS(2*(j-1)+1, n*(N+1)+m*(j-1)+4) = 1;
        ALOS(2*j, n*(j-1)+1) = -L3;
        ALOS(2*j, n*(j-1)+2) = -L4;
        ALOS(2*j, n*(N+1)+m*(j-1)+5) = -1;    
end
% Vectors for LOS cone constraints
bLOS = zeros(2*N, 1);
for i = 1:N
    bLOS(2*(i-1)+1) = 1;
    bLOS(2*(i-1)+2) = -1;
end

Ainit = zeros(n, numvar); Ainit(1:n, 1:n) = eye(n);
binit = x0vec;
% Concatenate all constraint matrices
Acon = [Aterm; ALOS];
bcon = [cb; bLOS];
Aeq = [Ainit; Adyn];
beq = [binit; zeros(size(Adyn,1),1)];
LB = llims;
UB = ulims;
f = zeros(numvar, 1);

H = zeros(numvar,numvar);
for i = 1:N
    H(n*(i-1)+1:n*i, n*(i-1)+1:n*i) = Q;
    H(n*(N+1)+1+m*(i-1):n*(N+1)+m*i, n*(N+1)+1+m*(i-1):n*(N+1)+m*i) = R;
end
H(n*N+1:n*(N+1), n*N+1:n*(N+1)) = P;


X = quadprog(H,f,Acon,bcon,Aeq,beq, LB, UB)
%X = quadprog(H,f,Acon,bcon,Aeq,beq)
LOSbl = -10^5.*ones(size(bcon));
Aconqpopt = [Aterm; ALOS; Ainit; Adyn];
uconqpopt = [ulims; bcon; beq];
lconqpopt = [llims; LOSbl; beq];
numvar = n*(N+1)+m*N
nclin = size(Aconqpopt,1)
ldA = nclin
ldH = size(H, 1)
leniw = 2*n+3;
lenw = 2*n^2+8*n+5*ldA;
x = zeros(numvar,1);
%writetestdat(fname_dat, ...
  %  Aconqpopt, H, lconqpopt, uconqpopt, cvec, x, x0vec)

%writeParamtest(fname_h,...
%    numvar, nclin, ldA, ldH, leniw, lenw)
 
 
 
 