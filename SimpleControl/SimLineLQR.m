clear all; close all; clc;
addpath(genpath('~/Software/tbxmanager'))
warning('off','MATLAB:nargchk:deprecated')

A = zeros(6,6);
A(1,3) = 1; A(2,4) = 1; A(5,6) = 1;
B = zeros(6,3);
B(3,1) = 1; B(4,2) = 1; B(6,3) = 1;
C = eye(6);
D = zeros(6,3);
sys = ss(A,B,C,D);

Q = 100.*eye(6); Q(3,3) = 10; Q(4,4) = 10;
R = 100.*eye(3); R(3,3) = 10;

[K,S,E] = lqr(sys,Q,R);

Pstruct.vxd = .03;
Pstruct.vyd = .03;
thetad = pi/2;
x0 = 0;
y0 = 0;
theta0 = pi/2;
y0 = [x0, y0, 0, 0, 0, 0, x0, y0, theta0]';
Pstruct.tdis0 = 2;
Pstruct.tdisf = 10;
Pstruct.udis = .01;
Pstruct.tdis = .01;
tspan = [0,10];
Aaug = zeros(9,9); 
Aaug(1:6,1:6) = A;
Aaug(7,1) = 1;
Aaug(8,2) = 1;
Aaug(9,5) = 1;
Baug = zeros(9,3);
Baug(1:6,:) = B;
Pstruct.A = Aaug;
Pstruct.B = Baug;
Pstruct.K = K;
Pstruct.umax = .4/18;
Pstruct.tmax = .78;
Pstruct.R = zeros(9,6);
Pstruct.R(7,1) = -1;
Pstruct.R(8,2) = -1;
Pstruct.R(9,5) = -1;
Pstruct.Ki = [.3 0 0; 0 .3 0; 0 0 .8];
[tout,yout] = ode45(@(t,x) spacecraft_dynfull(t,x,Pstruct),tspan,y0);

figure
plot(yout(:,1),yout(:,2))
hold all
plot(0:3,0:3)
xlim([0,max(yout(:,1))])
ylim([0,max(yout(:,2))])
figure
hold all
plot(tout,yout(:,1)-tout.*Pstruct.vxd)
plot(tout,yout(:,2)-tout.*Pstruct.vyd)
plot(tout,yout(:,5))
legend('x error','y error','theta error')
figure
hold all
plot(tout,yout(:,3))
plot(tout,yout(:,4))
plot(tout,yout(:,6))
legend('vx','vy','omega')

figure
hold all
plot(tout,yout(:,7))
plot(tout,yout(:,8))
plot(tout,yout(:,9))
legend('x int', 'y int','theta int')

Ktot = zeros(3,9);
Ktot(:,1:6) = K;
Ktot(:,7:9) = Pstruct.Ki;

folder_name = 'Line_filegen/';
fname_dat = 'Line_data.dat';
writeLineDat([folder_name fname_dat], Aaug, Baug, Ktot, Pstruct.R)
