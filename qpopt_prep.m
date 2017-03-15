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
u3_min = -20*pi/180;    % [rad/s] Desired rotation rate
u3_max = 20*pi/180;

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
    Adyn(n*(t-1)+1:n*(t-1)+n, n*t+1:n*t+n) = eye(n);
    Adyn(n*(t-1)+1:n*(t-1)+n, n*(N+1)+m*(t-1)+1:n*(N+1)+m*t) = Bd;
end

% Constructing vector of all dynamics constraints
bdyn = zeros(rowdim,1);
ldyn = bdyn; udyn = bdyn;

% Constructing vector of all state and thrust constraints
lxlims = zeros(numvar,1);
uxlims = zeros(numvar,1);
lulims = zeros(numvar,1);
uulims = zeros(numvar,1);

h_u = [u1_max u2_max u3_max 10^5 10^5 -u1_min -u2_min -u3_min -10^5 -10^5]';
h_x = [x1dot_max  x2dot_max  x1_max  x2_max  thetadot_max  theta_max...
      -x1dot_min -x2dot_min -x1_min -x2_min -thetadot_min -theta_min]';
  
for i=1:N+1
    luxlims(n*(i-1)+1:n*i,1) = h_x(1:n);
    lxlims(n*(i-1)+1:n*i,1) = h_x(n+1:2*n);
    if i<=N
        uulims(m*(i-1)+1:m*i,1) = h_u(1:m);
        lulims(m*(i-1)+1:m*i,1) = h_u(m+1:2*m);
    end
end
llims = [lxlims; lulims];
ulims = [uxlims; uulims];

% Constructing matrix of terminal constraints
nTc = length(cb);
Aterm = zeros(nTc, coldim);
Aterm(:, n*(N-1)+1:n*N) = cA;

% Constructing vector of terminal constraints
uterm = cb;
lterm = -10^6.*ones(size(uterm));

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
        ALOS(2*j, n*(j-1)+1) = L3;
        ALOS(2*j, n*(j-1)+1) = L4;
        ALOS(2*j, n*(N+1)+m*(j-1)+5) = 1;    
end
% Vectors for LOS cone constraints
lLOS = zeros(2*N, 1);
uLOS = zeros(2*N, 1);
for i = 1:N
    lLOS(2*(i-1)+1) = -10^5;
    lLOS(2*(i-1)+2) = 1;
    uLOS(2*(i-1)+1) = 1;
    uLOS(2*(i-1)+2) = 10^5;
end

% Concatenate all constraint matrices
Acon = [Adyn; Aterm; ALOS];
lcon = [llims; ldyn; lterm; lLOS];
ucon = [ulims; udyn; uterm; lLOS];
[nclin, numvar] = size(Acon)
ldA = ncong

H = zeros(numvar,numvar);
for i = 1:N
    H(n*(i-1)+1:n*i, n*(i-1)+1:n*i) = Q;
    H(n*(N+1)+1+m*(i-1):n*(N+1)+m*i, n*(N+1)+1+m*(i-1):n*(N+1)+m*i) = R;
end
H(n*N+1:n*(N+1), n*N+1:n*(N+1)) = P;

ldH = size(H,1)
leniw = 2*n+3
lenw = 2*n^2+8*n+5*ldA
cvec = zeros(ldH, 1);



