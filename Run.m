x0 =  3; y0 = 3; theta0 = -pi/4; vx0 = -.2; vy0 = .2; thetadot0 = -20*pi/180;
init = [x0 y0 theta0 vx0 vy0 thetadot0];
phi = pi/4; omega = 10*pi/180;
params.phi = phi;
params.omega = omega;
params.Ro = 550*10^3; % [m] orbital radius of target
params.mu = 3.986004418*10^14; % [m^3/s^2] standard gravitational parameter of earth
params.rp = .2; params.rtol = .2+.01;  params.rs = .2; params.gamma = pi/20;
params.Umax = 1; 
params.Tmax = 1;
params.Ts = .2;
params.N = 30; params.Nc = 15;
params.Qval = 10^3; params.Rval = 10^4; params.slackweight = 10^5;
params.vD = .7;
% slack_intercept and slack_slope are the intercept and slope respectively
% for determining the slack variable weight as a function of distance as
% the equation s = 10^(slack_slope*distance from goal point + slack_intercept)
% weight as a function of distance 
params.slack_intercept = 8; params.slack_slope = -0.75; 
params.slack_variable = 1; % If 1, slack weight will change with distance
params.eta = 1; params.betaHIGH = 1.5; params.betaLOW = 0.2;
[xtot, xtotI, utot, cost, time, ytot] = MPC_nonI_invset(init,params);
%%
x0 =  3; y0 = 3; theta0 = -pi/4; vx0 = .5; vy0 = .5; thetadot0 = -20*pi/180;
init = [x0 y0 theta0 vx0 vy0 thetadot0];
phi = pi/4; omega = 5*pi/180;
params.phi = phi;
params.omega = omega;
params.Ro = 850*10^3; % [m] orbital radius of target
params.mu = 3.986004418*10^14; % [m^3/s^2] standard gravitational parameter of earth
params.rp = .2; params.rtol = .2+.01;  params.rs = .2; params.gamma = pi/20;
params.Umax = .2;
params.Tmax = 1;
params.Ts = .2;
params.N = 30; params.Nc = 15;
params.Qval = 10^3; params.Rval = 10^4; params.slackweight = 10^5;
% slack_intercept and slack_slope are the intercept and slope respectively
% for determining the slack variable weight as a function of distance as
% the equation s = 10^(slack_slope*distance from goal point + slack_intercept)
% weight as a function of distance 
params.slack_intercept = 8; params.slack_slope = -0.75; 
params.slack_variable = 1; % If 1, slack weight will change with distance
params.eta = 1; params.betaHIGH = 1.5; params.betaLOW = 0.2;
[xtot, xtotI, utot, cost, time, ytot, slack] = MPC_Rotate_nonI(init,params);

%%
x0 =  3; y0 = 3; theta0 = -pi/4; vx0 = 0; vy0 = 0; thetadot0 = -20*pi/180;
init = [x0 y0 theta0 vx0 vy0 thetadot0];
phi = pi/4; omega = 0.6*pi/180;
params.phi = phi;
params.Ro = 450*10^3; % [m] orbital radius of target
params.mu = 3.986004418*10^14; % [m^3/s^2] standard gravitational parameter of earth
params.rp = .2; params.rtol = .2+.01;  params.rs = .2; params.gamma = pi/20;
params.Umax = .2;
params.Tmax = 1;
params.Ts = .2;
params.N = 15; params.Nc = 5;
params.Qval = 10^3; params.Rval = 10^3; params.slackweight = 10^5;
%params.Qval = 10^5; params.Rval = 10^2; params.slackweight = 10^8;
params.eta = .5; params.betaHIGH = 1.5; params.betaLOW = 0.2;
%params.betaHIGH = 100; params.betaLOW = -100; 
[xtot1, utot1, cost1, time1, ytot1] = MPC_Rotate_slack(init,params,phi,omega);

%%
x0 =  10; y0 = 0; z0 = 10; vx0 = 0; vy0 = 1; vz0 = 1; 
pred = 0; CWH = 1;
init = [x0 y0 z0 vx0 vy0 vz0];
params.nu = 0; params.phi = -pi/4;
params.Ro = 450*10^3; % [m] orbital radius of target
params.mu = 3.986004418*10^14; % [m^3/s^2] standard gravitational parameter of earth
params.rp = .2; params.rtol = .2+.01;  params.rs = .2; params.gamma = pi/20;
params.Umax = 3;
params.Ts = .2;
params.N = 20; params.Nc = 10;
[xtot, utot,cost] = MPC_3D_slack(init,params,pred,CWH);
%%
x0 = 5; y0 = 2; theta0 =-.5; vx0 = 0; vy0 = 1; omega0 = 0; phi = .5; omegaD = 10*pi/180;
init = [x0 y0 theta0 vx0 vy0 omega0];
params.rp = .2; params.rtol = .2+.01;  params.rs = .2; params.gamma = pi/20;
params.Ts = .2;
params.eta = 1; params.etaR = 1; params.beta1 = .05; params.beta2 = 2; 
params.beta3 = .5; params.beta4 = omegaD*1.7; params.beta5 = -omegaD*1.2;
params.N = 20; params.Nc = 10;
thetaD = phi; 
x0vec = [x0; y0; theta0-phi; vx0; vy0; omega0];
[xtot, utot,ytot] = MPC_Quad(init,params,phi);
xD = [0 0 phi 0 0 0];
%%
x0 = 5; y0 = 2; theta0 =-.5; vx0 = 0; vy0 = 1; omega0 = 0; phi = .5; omegaD = 10*pi/180;
init = [x0 y0 theta0 vx0 vy0 omega0];
params.rp = .2; params.rtol = .2+.05;  params.rs = .2; params.gamma = 10*pi/180;
params.Ts = .2;
params.eta = 1; params.etaR = 1; params.beta1 = .05; params.beta2 = 2; 
params.beta3 = .5; params.beta4 = omegaD*1.7; params.beta5 = -omegaD*1.2;
params.N = 20; params.Nc = 10;
thetaD = phi; 
x0vec = [x0; y0; theta0-phi; vx0; vy0; omega0];
[xtot, utot,ytot] = CVX_coneandALLvel_editBOUNDS(init,params,phi);
xD = [0 0 phi 0 0 0];

%%
x0 = 5; y0 = 2; theta0 =-.5; vx0 = 0; vy0 = 1; omega0 = 0; phi = .5; omegaD = 20*pi/180;
init = [x0 y0 theta0 vx0 vy0 omega0];
params.rp = .2; params.rtol = .2+.05; params.rs = .2; params.gamma = 10*pi/180; 
params.Ts = .2; 
params.eta = 1; params.etaR = .2; params.beta1 = .05; params.beta2 = 1.5; 
params.beta3 = .2; params.beta4 = omegaD*1.2; params.beta5 = -omegaD*.8; 
params.rho = 10^10;
params.N = 20; params.Nc = 10;

eta = 1; beta1 = .05; beta2 = 2; beta3 = .5; beta4 = omegaD*3.5; rho = 10^10;
beta5 = -omegaD*.5; 
beta6 = phi*1.8; beta7 = -phi/10;
thetaD = phi; 
%phi = .5;
%x0 = 5; y0 = 2; 
%theta0 = 90*pi/180; vx0 = 0; vy0 = 0; omega0 = 0; 
x0vec = [x0; y0; theta0-phi; vx0; vy0; omega0];
[xtot, utot,ytot] = CVX_coneandALLvel_changeQ(x0,y0,theta0,vx0,vy0,omega0,phi,omegaD);
xD = [0 0 phi 0 0 0];

%%
animate_path(1, xtot, xD, rp, rs, rtol, gamma, 1,0,'highspeed_allVELcon_rotation.avi')
%%
animate_path_nooffset(xtot,rp,rs,rtol,gamma,1,phi)
%%
x0 = 5; y0 = 2; theta0 =-.5; vx0 = 0; vy0 = 1; omega0 = 0; phi = .5; 
omegaDmat = linspace(-20*pi/180,20*pi/180,10);
for i=1:10
    omegaD = omegaDmat(i);
    if omegaD<0
        theta
x0vec = [x0; y0; theta0-phi; vx0; vy0; omega0];
[xtot, utot,ytot] = CVX_coneandALLvel_WORKING(x0,y0,theta0,vx0,vy0,omega0,phi,omegaD);

end





