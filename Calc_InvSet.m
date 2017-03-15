function [cA, cb] = Calc_InvSet(params, A, B)

x1dot_min = -0.5;       % m/s (true limit is ~0.3)
x1dot_max = 0.5; 
x2dot_min = -0.5;
x2dot_max = 0.5;
% x1_min = 0.1;         % meters (true limits)
% x1_max = 3.5;   
% x2_min = 0.1;
% x2_max = 2.5;
x1_min = -10;           % meters (relaxed limits)
x1_max = 10;   
x2_min = -10;
x2_max = 10;
thetadot_min = -pi;     % rad/s (made up)
thetadot_max = pi; 
theta_min = -50*pi;     % rad (in principle, should have no limits)
theta_max = 50*pi;
u1_min = -.2;            % Two thrusters, normalized
u1_max = .2;             % Cut actuation limit in half due to 50% actuation duty cycle
u2_min = -.2;
u2_max = .2;
u3_min = -.5;    % [rad/s] Desired rotation rate
u3_max = .5;

% Set final state constraint to maximally invariant set
model = LTISystem('A',A,'B',B);
model.x.min = [x1_min x2_min theta_min x1dot_min x2dot_min thetadot_min]';
model.x.max = [x1_max x2_max theta_max x1dot_max x2dot_max thetadot_max]';
model.u.min = [u1_min u2_min u3_min];
model.u.max = [u1_max u2_max u3_max];
invset = model.invariantSet('maxIterations',40);
cA = invset.A; cb = invset.b;