% Generates the points on the spherical target that the cosntraint plane
% must intersect with

theta = pi/4; phi = -pi/4; rho = 5*pi/180;
rp = 0.2; % [m] radius of target
R1 = [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];
R2 = [cos(phi) 0 -sin(phi); 0 1 0; sin(phi) 0 cos(phi)];
R = R2*R1;
Rinv = inv(R);
[X,Y,Z] = sphere;

rvec = Rinv*[rp; 0; 0];
Xt = X.*rp; Yt = Y.*rp; Zt = Z.*rp;
Xp = X.*rp/40; Yp = Y.*rp/40; Zp = Z.*rp/40;
surf(Xt,Yt,Zt,'Facecolor', [0 1 0])
hold all
surf(Xp+rvec(1), Yp+rvec(2), Zp+rvec(3),'Facecolor',[1 0 0])
xlabel('X axis')
ylabel('Y axis')
zlabel('Z axis')

view([120 20])

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

p1avec = Rp1ainv*[rp; 0; 0];
p1bvec = Rp1binv*[rp; 0; 0];
p2avec = Rp2ainv*[rp; 0; 0];
p2bvec = Rp2binv*[rp; 0; 0];
p3avec = Rp3ainv*[rp; 0; 0];
p3bvec = Rp3binv*[rp; 0; 0];
p4avec = Rp4ainv*[rp; 0; 0];
p4bvec = Rp4binv*[rp; 0; 0];

surf(Xp+p1avec(1), Yp+p1avec(2), Zp+p1avec(3),'Facecolor',[1 0 0])
surf(Xp+p1bvec(1), Yp+p1bvec(2), Zp+p1bvec(3),'Facecolor',[1 0 0])
surf(Xp+p2avec(1), Yp+p2avec(2), Zp+p2avec(3),'Facecolor',[1 0 0])
surf(Xp+p2bvec(1), Yp+p2bvec(2), Zp+p2bvec(3),'Facecolor',[1 0 0])
surf(Xp+p3avec(1), Yp+p3avec(2), Zp+p3avec(3),'Facecolor',[1 0 0])
surf(Xp+p3bvec(1), Yp+p3bvec(2), Zp+p3bvec(3),'Facecolor',[1 0 0])
surf(Xp+p4avec(1), Yp+p4avec(2), Zp+p4avec(3),'Facecolor',[1 0 0])
surf(Xp+p4bvec(1), Yp+p4bvec(2), Zp+p4bvec(3),'Facecolor',[1 0 0])


