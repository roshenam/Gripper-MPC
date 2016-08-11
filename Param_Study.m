% Run simulations across variations of Q and R and N with varying x, y, xdot,
% ydot to check region of attraction, cost variations, and final conditions
phi = pi/4; omega = 2*pi/180;
params.phi = pi/4; % Phi is the angle of the port from the defined horizontal
params.omega = 2*pi/180; % Omega is the rotation in rad/sec of the target
params.rp = .2; params.rtol = .2+.01;  params.rs = .2; params.gamma = 20*pi/180;
% rp is the radius of the target, rs is the radius of the spacecraft, gamma
% is the half angle of the constraint cone
params.Umax = .5;
% Max thrust in Newtons/kg of the free flyer
params.Tmax = 1;
% Max torque in N*m/(kg*m^2)
params.Ts = .2;
% Discretization constant (0.2 seconds)
params.Qval = 10^4; params.Rval = 10^4;
% Weights on states, control inputs, and slack variables
params.betaHIGH = -1.5;
params.eta = 1;
params.tantol = .4;
params.N = 30;
params.Nsim = 50;
states = [3 .2 0 0; 3 .2 -.5 -.5; 3 .2 -.5 .5; 3 .2 .5 -.5; 3 .2 .5 .5];
% Constants defining the shape of the bounding function for normal
% component of hte velocity
Xtot = cell(4,1); Utot = cell(4,1); XtotnonI = cell(4,1); 
UtotnonI = cell(4,1); Final = cell(4,1);
Init = cell(4,1);

for i=1:4
x0 =  states(i,1); y0 = states(i,2); theta0 = -pi/4; vx0 = states(i,3); vy0 = states(i,4); thetadot0 = 5*pi/180;
Rmat = [cos(phi) -sin(phi); sin(phi) cos(phi)];
rI = Rmat*[x0;y0];
vI = Rmat*[vx0;vy0];
init = [rI(1) rI(2) theta0 vI(1) vI(2) thetadot0];
x0vec = [x0; y0; nu0; vx0; vy0; nu0dot];
[xtot, utot, xtotnonI, utotnonI, final] = Inv_Set_Sim(x0vec, params);
Xtot{i} = xtot; Utot{i} = utot; XtotnonI{i} = xtotnonI; 
UtotnonI{i} = utotnonI; Final{i} = final;
Init{i} = init;
end

%%
for i = 1:4
    Animate(Init{i},params,Xtot{i},Xtot{i}(end,:),0,'no',1,1)
    pause
end
