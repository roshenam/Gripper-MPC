% x, y, xdot, ydot, t, tdot
omega = 10*pi/180;
Ts = .2;
Umax = 1;
A = zeros(6,6); A(1,3) = 1; A(2,4) = 1; A(5,6) = 1;
A(3,1) = omega^2; A(3,4) = 2*omega; 
A(4,2) = omega^2; A(4,3) = -2*omega;
B = zeros(6,2); B(3,1) = 1; B(4,2) = 1;
C = zeros(6,6); D = zeros(6,2);
sysD = ss(A,B,C,D);
sysD = c2d(sysD,Ts);
vD = .5;
gamma = 10*pi/180;
x0vec = [6; .5; .2+vD; .2; 0; 1];
Q = 1000.*eye(6); Q(5,5) = 0; Q(6,6) = 0;
R = 100.*eye(2);
[~,P,~] = lqrd(A(1:4,1:4),B(1:4,:),Q(1:4,1:4),R,Ts);
Pbig = zeros(6,6);
Pbig(1:4,1:4) = P;
Abig = sysD.a; Bbig = sysD.b; 
Cbig = zeros(2,6); Cbig(1,1) = -tan(gamma); Cbig(1,2) = 1; Cbig(1,5) = tan(gamma)*vD;
Cbig(2,1) = -tan(gamma); Cbig(2,2) = -1; Cbig(2,5) = tan(gamma)*vD;
Dbig = zeros(2,2);
Ymax = [0;0];
rp = .5; rs = .5;
counter = 0;
time = [];
utot = [];
xtot = x0vec;
cost = [];
n = 6; m = 2; p = 2; N = 15;
%while norm([x0vec(1)-vD*counter*Ts,x0vec(2)])>=(rp+rs)
for j=1:30
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
    %max(U(3,:)') <= Tmax;
    %min(U(3,:)') >= -Tmax;
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
    %ytot = [ytot Ymax];
    currcost = cvx_optval;
    cost = [cost; currcost];
    if counter>200
        disp(['More than 200 iterations. Stopping simulation.'])
        break
    end
    
end
disp('Simulation Complete')
xreal = xtot(1,:) - vD.*xtot(5,:);
vreal = xtot(2,:) - vD;
