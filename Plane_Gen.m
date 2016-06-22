%% Generates points on the target to intersect with constraint planes

theta = -pi/4; phi = -pi/4; rho = 10*pi/180;
rp = 0.2; % [m] radius of target
R1 = [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];
R2 = [cos(phi) 0 -sin(phi); 0 1 0; sin(phi) 0 cos(phi)];
R = R2*R1;
Rinv = inv(R);
[X,Y,Z] = sphere;

rvec = Rinv*[rp; 0; 0];
Xt = X.*rp; Yt = Y.*rp; Zt = Z.*rp;
Xp = X.*rp/40; Yp = Y.*rp/40; Zp = Z.*rp/40;

% Rotations to generate points at NSWE corners of circle
Ry = [cos(rho) 0 -sin(rho); 0 1 0; sin(rho) 0 cos(rho)];
Ryn = [cos(-rho) 0 -sin(-rho); 0 1 0; sin(-rho) 0 cos(-rho)];
Rz = [cos(rho) sin(rho) 0; -sin(rho) cos(rho) 0; 0 0 1];
Rzn = [cos(-rho) sin(-rho) 0; -sin(-rho) cos(-rho) 0; 0 0 1];
Rp1a = Ry*R; Rp1ainv = inv(Rp1a); 
Rp1b = Ryn*R; Rp1binv = inv(Rp1b); 
Rp2a = Rz*R; Rp2ainv = inv(Rp2a); 
Rp2b = Rzn*R; Rp2binv = inv(Rp2b); 

% Rotations to generate points at intermediary corners of circle
Rx = [1 0 0; 0 cos(pi/4) sin(pi/4); 0 -sin(pi/4) cos(pi/4)];
Rxn = [1 0 0; 0 cos(-pi/4) sin(-pi/4); 0 -sin(-pi/4) cos(-pi/4)];
Rs = [cos(rho) 0 -sin(rho); 0 1 0; sin(rho) 0 cos(rho)];
Rn = [cos(-rho) 0 -sin(-rho); 0 1 0; sin(-rho) 0 cos(-rho)];

Rp3a = Rs*Rx*R; Rp3ainv = inv(Rp3a);
Rp3b = Rn*Rx*R; Rp3binv = inv(Rp3b);
Rp4a = Rs*Rxn*R; Rp4ainv = inv(Rp4a);
Rp4b = Rn*Rxn*R; Rp4binv = inv(Rp4b);

pS = Rp1ainv*[rp; 0; 0];
pN = Rp1binv*[rp; 0; 0];
pE = Rp2ainv*[rp; 0; 0];
pW = Rp2binv*[rp; 0; 0];
pSE = Rp3ainv*[rp; 0; 0];
pNW = Rp3binv*[rp; 0; 0];
pSW = Rp4ainv*[rp; 0; 0];
pNE = Rp4binv*[rp; 0; 0];

%% Plot points to make sure the code is making sense
target = surf(Xt,Yt,Zt,'Facecolor', [0 1 0]);
hold all
surf(Xp+rvec(1), Yp+rvec(2), Zp+rvec(3),'Facecolor',[1 0 0])
% surf(Xp+pS(1), Yp+pS(2), Zp+pS(3),'Facecolor',[1 0 0])
% surf(Xp+pN(1), Yp+pN(2), Zp+pN(3),'Facecolor',[1 0 0])
% surf(Xp+pE(1), Yp+pE(2), Zp+pE(3),'Facecolor',[1 0 0])
% surf(Xp+pW(1), Yp+pW(2), Zp+pW(3),'Facecolor',[1 0 0])
% surf(Xp+pSE(1), Yp+pSE(2), Zp+pSE(3),'Facecolor',[1 0 0])
% surf(Xp+pNW(1), Yp+pNW(2), Zp+pNW(3),'Facecolor',[1 0 0])
% surf(Xp+pSW(1), Yp+pSW(2), Zp+pSW(3),'Facecolor',[1 0 0])
% surf(Xp+pNE(1), Yp+pNE(2), Zp+pNE(3),'Facecolor',[1 0 0])

xlabel('X axis')
ylabel('Y axis')
zlabel('Z axis')

view([90+theta*180/(2*pi) -phi*180/(2*pi)])

%% Plot 8 planes to approximate the cone
X = [0 0 0 0 0 0 0 0; 2*pN(1) 2*pNE(1) 2*pE(1) 2*pSE(1) 2*pS(1) 2*pSW(1) 2*pW(1) 2*pNW(1);
    2*pNE(1) 2*pE(1) 2*pSE(1) 2*pS(1) 2*pSW(1) 2*pW(1) 2*pNW(1) 2*pN(1)];
Y = [0 0 0 0 0 0 0 0; 2*pN(2) 2*pNE(2) 2*pE(2) 2*pSE(2) 2*pS(2) 2*pSW(2) 2*pW(2) 2*pNW(2);
    2*pNE(2) 2*pE(2) 2*pSE(2) 2*pS(2) 2*pSW(2) 2*pW(2) 2*pNW(2) 2*pN(2)];
Z = [0 0 0 0 0 0 0 0; 2*pN(3) 2*pNE(3) 2*pE(3) 2*pSE(3) 2*pS(3) 2*pSW(3) 2*pW(3) 2*pNW(3);
    2*pNE(3) 2*pE(3) 2*pSE(3) 2*pS(3) 2*pSW(3) 2*pW(3) 2*pNW(3) 2*pN(3)];
C = [0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5;
     1 1 1 1 1 1 1 1;
     0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3];
h = fill3(X,Y,Z,C);
alpha(h, .5)

%% Compute normal vectors for each plane
p1 = cross(pN, pNE); p2 = cross(pNE, pE); 
p3 = cross(pE, pSE); p4 = cross(pSE, pS);
p5 = cross(pS, pSW); p6 = cross(pSW, pW); 
p7 = cross(pW, pNW); p8 = cross(pNW, pN);

