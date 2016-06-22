function [xtot, utot] = MPC_3D(init,params,phi,nu,pred)
% If pred is zero, there will be no prediction of the LOS cone constraints
% in the model. If pre is one, there will be prediction of the LOS cone
% constraints. Rotation for docking port is first rotate about the z axis
% with angle nu, then rotate about the new x axis with angle phi.

% 6/20 First doing without any attitude control to avoid having to
% linearize equations. Translation and rotation separate.

% Extracting initial conditions and parameters from inputs
x0 = init(1); y0 = init(2); z0 = init(3);  
vx0 = init(4); vy0 = init(5); vz0 = init(6);
rp = params.rp; rs = params.rs; gamma = params.gamma;
Ts = params.Ts;  
N = params.N; Nc = params.Nc; 

x0vec = [x0; y0; z0; vx0; vy0; vz0];

% Creating dynamic system 
A = zeros(6,6); A(1,4) = 1; A(2,5) = 1; A(3,6) = 1; 
B = zeros(6,4); B(4,1) = 1; B(5,2) = 1; B(6,3) = 1;
C = eye(6); D = zeros(6,4);

% Upper bounds on thrust inputs and slack variables
Umax = [0.2 0.2 0.2 10^10]';
Umin = -Umax;
sysD = ss(A,B,C,D);
sysD = c2d(sysD,Ts);
% Creating weighting matrices for cost function
Q = 10^3.*eye(6); R = 100.*eye(4); 
R(4,4) = 10^6;
% Solving infinite horizon unconstrained LQR for stability enforcement
[K,P,~] = lqrd(A,B(:,1:3),Q,R(1:3,1:3),Ts);

n = 6; m = 4;
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
