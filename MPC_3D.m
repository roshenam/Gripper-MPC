function [xtot, utot] = MPC_3D(init,params,pred,CWH)
% If pred is zero, there will be no prediction of the LOS cone constraints
% in the model. If pre is one, there will be prediction of the LOS cone
% constraints. Rotation for docking port is first rotate about the z axis
% with angle nu, then rotate about the new x axis with angle phi. If CWH is
% zero, there will be no use of CWH dynamics. If CWH is 1, CWH dynamics
% will be taken into account and orbital parameters must be specified in
% the params structure. Constraints are implemented as hard. MPC_3D_slack
% and MPC_3D_cost for soft implementations of constraints

% 6/20 First doing without any attitude control to avoid having to
% linearize equations. Translation and rotation separate.
% 6/23 Finished CWH dynamics implementation.

% Extracting initial conditions and parameters from inputs
x0 = init(1); y0 = init(2); z0 = init(3);
vx0 = init(4); vy0 = init(5); vz0 = init(6);
phi = params.phi; nu = params.nu;
rp = params.rp; rs = params.rs;
Ts = params.Ts;
N = params.N; Nc = params.Nc;

x0vec = [x0; y0; z0; vx0; vy0; vz0];

% Creating dynamic system
if CWH
    % Accounting for CWH dynamics
    mu = params.mu; Ro = params.Ro;
    n = sqrt(mu/(Ro^3));
    A = zeros(6,6); A(1:3,4:6) = eye(3);
    A(4,1) = 3*n^2; A(4,5) = 2*n;
    A(5,4) = -2*n; A(6,3) = -n^2;
else
    % Not accounting for CWH dynamics
    A = zeros(6,6); A(1,4) = 1; A(2,5) = 1; A(3,6) = 1;
end
B = zeros(6,3); B(4,1) = 1; B(5,2) = 1; B(6,3) = 1;

C = zeros(8,6); D = zeros(8,3);

% C is a matrix of the weights on the constrained outputs due to the 8
% planes approximating the cone
Cs = Plane_Gen(init, params,phi,nu);
for k=1:8
    C(k,1) = Cs(1,k);
    C(k,2) = Cs(2,k);
    C(k,3) = Cs(3,k);
end

% Upper bounds on thrust inputs and slack variables
Umax = [6.8 1 1]';
Umin = -Umax;
sysD = ss(A,B,C,D);
sysD = c2d(sysD,Ts);


% Creating weighting matrices for cost function
Q = 10^3.*eye(6); R = 100.*eye(3);
% Solving infinite horizon unconstrained LQR for stability enforcement
[K,P,~] = lqrd(A,B(:,1:3),Q,R(1:3,1:3),Ts);

n = 6; m = 3;
utot = [];
ytot = [];
xtot = x0vec;
counter = 0;
while norm(x0vec(1:2))>=(rp+rs)
    counter = counter + 1;
    disp(['Running optimization ',num2str(counter),', distance from origin is ',num2str(norm(x0vec(1:2)))])
    cvx_begin quiet
    variables X(n,N+1) U(m,N)
    X(:,2:N+1) == sysD.a*X(:,1:N) + sysD.b*U;
    X(:,1) == x0vec;
    max(C*X(:,1:Nc)) <= 0;
    max(U(1,:).^2 + U(2,:).^2 + U(3,:).^2) <= Umax(1)^2;
    
    % CVX returns that a portion of this is illegal ({convex}^.5), although
    % the whole thing is convex (cone shape) so not sure how to get around
    % this. For now, ditch this and go with linear approximations of the
    % cone like in the paper.
    %     max(pow_p((X(3,:).*sin(phi) + X(2,:).*cos(nu)*cos(phi) - X(1,:).*cos(phi)*sin(nu)).^2 +...
    %     (X(1,:).*cos(nu) + X(2,:).*sin(nu)).^2,0.5) - tan(beta).*(X(3,:).*cos(phi) -...
    %     X(2,:).*cos(nu)*sin(phi) + X(1,:).*sin(nu)*sin(phi)) + U(4,:)) <= 0;
    
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
    
end
disp('Simulation Complete')
