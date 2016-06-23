function [] = Animate_3D(init, params, phi, nu, x, animate,filename)
if animate
v = VideoWriter(filename);
v.FrameRate = 6;
open(v);
end
gamma = params.gamma; rp = params.rp; 
rp = 0.2; % [m] radius of target
R1 = [cos(nu) sin(nu) 0; -sin(nu) cos(nu) 0; 0 0 1];
R2 = [cos(phi) 0 -sin(phi); 0 1 0; sin(phi) 0 cos(phi)];
R = R2*R1;
Rinv = inv(R);
[X,Y,Z] = sphere;

rvec = Rinv*[rp; 0; 0];
Xt = X.*rp; Yt = Y.*rp; Zt = Z.*rp;
Xp = X.*rp/40; Yp = Y.*rp/40; Zp = Z.*rp/40;

% Rotations to generate points at NSWE corners of circle
Ry = [cos(gamma) 0 -sin(gamma); 0 1 0; sin(gamma) 0 cos(gamma)];
Ryn = [cos(-gamma) 0 -sin(-gamma); 0 1 0; sin(-gamma) 0 cos(-gamma)];
Rz = [cos(gamma) sin(gamma) 0; -sin(gamma) cos(gamma) 0; 0 0 1];
Rzn = [cos(-gamma) sin(-gamma) 0; -sin(-gamma) cos(-gamma) 0; 0 0 1];
Rp1a = Ry*R; Rp1ainv = inv(Rp1a); 
Rp1b = Ryn*R; Rp1binv = inv(Rp1b); 
Rp2a = Rz*R; Rp2ainv = inv(Rp2a); 
Rp2b = Rzn*R; Rp2binv = inv(Rp2b); 

% Rotations to generate points at intermediary corners of circle
Rx = [1 0 0; 0 cos(pi/4) sin(pi/4); 0 -sin(pi/4) cos(pi/4)];
Rxn = [1 0 0; 0 cos(-pi/4) sin(-pi/4); 0 -sin(-pi/4) cos(-pi/4)];
Rs = [cos(gamma) 0 -sin(gamma); 0 1 0; sin(gamma) 0 cos(gamma)];
Rn = [cos(-gamma) 0 -sin(-gamma); 0 1 0; sin(-gamma) 0 cos(-gamma)];

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

%Plot points to make sure the code is making sense
target = surf(Xt,Yt,Zt,'Facecolor', [0 1 0]);
hold all
surf(Xp+rvec(1), Yp+rvec(2), Zp+rvec(3),'Facecolor',[1 0 0])
surf(Xp+pS(1), Yp+pS(2), Zp+pS(3),'Facecolor',[1 0 0])
surf(Xp+pN(1), Yp+pN(2), Zp+pN(3),'Facecolor',[1 0 0])
surf(Xp+pE(1), Yp+pE(2), Zp+pE(3),'Facecolor',[1 0 0])
surf(Xp+pW(1), Yp+pW(2), Zp+pW(3),'Facecolor',[1 0 0])
surf(Xp+pSE(1), Yp+pSE(2), Zp+pSE(3),'Facecolor',[1 0 0])
surf(Xp+pNW(1), Yp+pNW(2), Zp+pNW(3),'Facecolor',[1 0 0])
surf(Xp+pSW(1), Yp+pSW(2), Zp+pSW(3),'Facecolor',[1 0 0])
surf(Xp+pNE(1), Yp+pNE(2), Zp+pNE(3),'Facecolor',[1 0 0])
plot3(init(1), init(2), init(3),'ko')
xlabel('X axis')
ylabel('Y axis')
zlabel('Z axis')

%view([90+nu*180/(2*pi) -phi*180/(2*pi)])
view([51 17])

% Plot 8 planes to approximate the cone
X = 8.*[0 0 0 0 0 0 0 0; pN(1) pNE(1) pE(1) pSE(1) pS(1) pSW(1) pW(1) pNW(1);
    pNE(1) pE(1) pSE(1) pS(1) pSW(1) pW(1) pNW(1) pN(1)];
Y = 8.*[0 0 0 0 0 0 0 0; pN(2) pNE(2) pE(2) pSE(2) pS(2) pSW(2) pW(2) pNW(2);
    pNE(2) pE(2) pSE(2) pS(2) pSW(2) pW(2) pNW(2) pN(2)];
Z = 8.*[0 0 0 0 0 0 0 0; pN(3) pNE(3) pE(3) pSE(3) pS(3) pSW(3) pW(3) pNW(3);
    pNE(3) pE(3) pSE(3) pS(3) pSW(3) pW(3) pNW(3) pN(3)];
C = [0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5;
     1 1 1 1 1 1 1 1;
     0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3];
h = fill3(X,Y,Z,C);
alpha(h, .5);


freeflyer = plot3(init(1),init(2),init(3),'ro','markersize',3);
path = line(init(1),init(2),init(3),'color','b');

frame = getframe;
if animate
writeVideo(v,frame);
end
for i=2:length(x)
    set(freeflyer, 'XData', x(1,i), 'YData', x(2,i), 'ZData', x(3,i));
    set(path,'XData',x(1,1:i),'YData',x(2,1:i), 'ZData',x(3,1:i));
    drawnow  
    pause(.1)
    if animate
    frame = getframe;
    writeVideo(v,frame);
    end
end
if animate
close(v);
end