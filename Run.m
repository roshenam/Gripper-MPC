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





