% Generates the structure array of parameters that's an input for
% Gen_InvSet

% rt is radius of target, rtol should be radius of target plus a little bit
% to ensure feasibility, rs is radius of spacecraft aka free flyer, w is
% the distance from the edge of the spacecraft to the grasping surface of
% the gripper
params.rt = .15; params.rtol = params.rt+.02;  params.rs = .15;
params.w = .1;
% theta_c is the worst  case offset of the angle from the desired angle at
% the end, it's an estimate
params.theta_c = 2*pi/180;
% dmax is the worst case offset we can have. this should be a function of
% omega! we need to pull this from the experimental data. for now since w
% don't have the thrust for this to work at much more than 5-10 deg/sec, we
% can see from the plot in Matt's paper that a max offset of 5 cm should be
% ok
params.dmax = .05;
% rt is the radius of the target, rs is the radius of the spacecraft
params.Umax = .2;
% Max thrust in Newtons/kg of the free flyer
params.Tmax = .3; % 50% of stall torque
% Max torque in N*m/(kg*m^2)
params.Ts = 0.5;
% Discretization constant (in seconds)
params.Qval = 10^3; params.Rval = 10^2;
params.Qval2 = 10^2;
params.slackweight = 10^6;

% Do not consider omegas greater than 4 deg/sec with the thrust max values
% given above, the set will most likely be empty because we don't have
% enough thrust to be able to stay in the feasible region

% Example call to Gen_InvSet
[Acon, bcon, flag] = Gen_InvSet( params, 4*pi/180);
