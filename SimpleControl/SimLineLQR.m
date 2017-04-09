clear all; 
%close all; 
clc;
addpath(genpath('~/Software/tbxmanager'))
warning('off','MATLAB:nargchk:deprecated')

A = zeros(6,6);
A(1,3) = 1; A(2,4) = 1; A(5,6) = 1;
B = zeros(6,3);
B(3,1) = 1; B(4,2) = 1; B(6,3) = 1;
C = eye(6);
D = zeros(6,3);
sys = ss(A,B,C,D);

% Q = 100.*eye(6); Q(3,3) = 10; Q(4,4) = 10;
% R = 100.*eye(3); R(3,3) = 10;

Q = eye(6); Q(3,3) = 100; Q(4,4) = 100; Q(5,5) = 10; Q(6,6) = 1;
R = 10^2.*eye(3);

[K,S,E] = lqr(sys,Q,R);

Pstruct.vxd = .03;
Pstruct.vyd = .03;
thetad = pi/2;
x0 = .01;
y0 = 0;
theta0 = pi/2;
y0 = [x0, y0, 0, 0, theta0, 0, x0, y0, theta0]';
Pstruct.tdis0 = 2;
Pstruct.tdisf = 4;
Pstruct.udis = 0;
Pstruct.torquedis = 0;
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
Pstruct.umax = (.4/18.08)/sqrt(2);
Pstruct.tmax = 1.48;
Pstruct.R = zeros(9,6);
Pstruct.R(7,1) = -1;
Pstruct.R(8,2) = -1;
Pstruct.R(9,5) = -1;
Pstruct.Ki = [0 0 0; 0 0 0; 0 0 0];
[tout,yout] = ode45(@(t,x) spacecraft_dynfull(t,x,Pstruct),tspan,y0);

figure(1)
hold all
plot(yout(:,1),yout(:,2))
plot([0,max(yout(:,1))],[0,max(yout(:,2))])
hold all
figure(2)
hold all
plot(yout(:,5))
% plot(tout,yout(:,1)-tout.*Pstruct.vxd)
% plot(tout,yout(:,2)-tout.*Pstruct.vyd)
% plot(tout,yout(:,5))
% legend('x error','y error','theta error')
figure(3)
hold all
plot(tout,yout(:,3))
plot(tout,yout(:,4))
plot(tout,yout(:,6))
legend('vx','vy','omega')

% figure
% hold all
% plot(tout,yout(:,7))
% plot(tout,yout(:,8))
% plot(tout,yout(:,9))
% legend('x int', 'y int','theta int')

Ktot = zeros(3,9);
Ktot(:,1:6) = K;
Ktot(:,7:9) = Pstruct.Ki;
Ktot
Pstruct.R


% plot(yout(:,3))
% hold all
% plot(yout(:,4))
% figure
% plot(yout(:,1)-Pstruct.vxd.*tout)

%Ktot = zeros(3,8);
%Ktot(:,1:6) = K;
%Ktot(:,7:9) = Pstruct.Ki;


folder_name = 'Line_filegen/';
fname_dat = 'Line_data.dat';
writeLineDat([folder_name fname_dat], K)
