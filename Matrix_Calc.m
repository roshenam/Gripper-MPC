
%% Free-flyer specifications
in2m = 0.0254;      % Inch to meter conversion

% Physical parameters
mass = 16.3;        % Mass (kg)
F = 0.2;            % Thruster force (single thruster) (N)
l = 9.25/2*in2m;    % Thruster arm (m)
Mthrust = 4*F*l;    % Dual thruster torque (N*m)

% Guesswork
mass_w = 3.64;      % Mass of wheel
mass_ff = mass - mass_w;    % Mass of free-flyer
r_w = 5*in2m;       % Wheel radius
r_ff = 5.5*in2m;    % Free-flyer radius
Iw = 0.5*mass_w*r_w^2;
If = 0.5*mass_ff*r_ff^2;     % Approximate as cylinder

J=Iw+If;            % Inertia, full assembly (kg*m^2)
alpha = Mthrust/J;  % Rotational acceleration due to dual thruster torque

% Motor Parameters
Ncrev = 64;                                 % [counts/rev] Counts per motor revolution
Ngr = 18.75;                                % Gear ratio
Nc_absmax = 2047;                           % [counts] Absolute maximum number of counts measurable per PID period
rpms_max = 250;                             % [RPM] No load shaft RPM
fs_max = rpms_max/60;                       % [rev/s] No load shaft rotation frequency
ws_max = fs_max*2*pi;                       % [rad/s] No load shaft angular velocity
wm_max = ws_max*Ngr;                        % [rad/s] No load motor angular velocity
wm_mid = wm_max/2;                          % [rad/s] Motor angular velocity when free-flyer angular velocity is zero
ws_mid = ws_max/2;                          % [rad/s] Shaft angular velocity when free-flyer angular velocity is zero

ozin2Nm = 0.007062;                         % Convert ounce-inches to Newton-meters
tau_max = 84*ozin2Nm;                       % [N*m] Maximum torque delivered by motor
Ndc_max = 600;                              % Maximum duty cycle target

T_PID = 0.02;                               % [s] Motor controller PID Period
If = 0.0816;                                % [kg*m^2] Free-flyer moment of inertia (not including reaction wheel)
Iw = 0.0294;                                % [kg*m^2] Reaction wheel moment of inertia

Nc_max = Ncrev*Ngr*fs_max*T_PID             % [counts] Maximum number of counts measurable per PID period
max2absmax = Nc_absmax/Nc_max;

% If free-flyer is floating
ws_absmax = 2*pi*Nc_absmax/(Ncrev*Ngr*T_PID);   % [rad/s] Absolute maximum shaft angular velocity (assuming maximum is not reached)
wf_absmax = ws_absmax/(If/Iw-1);                % [rad/s] Absolute maximum free-flyer angular velocity
ff_absmax = wf_absmax/(2*pi);                   % [rev/s] Absolute maximum free-flyer rotation frequency

wf_max = ws_max/(If/Iw-1);                  % [rad/s] Maximum free-flyer angular velocity
ff_max = wf_max/(2*pi);                     % [rev/s] Maximum free-flyer rotation frequency

% If free-flyer is stationary
ww_absmax = ws_absmax;                      % [rad/s] Absolute maximum wheel angular velocity
fw_absmax = ww_absmax/(2*pi);               % [rev/s] Absolute maximum wheel rotation frequency
fw_max = fs_max;                            % [rev/s] Maximum rotation frequency

Kp = 4;                                             % Proportional PID gain
Ki = 0.01;                                          % Integral PID gain
Ky = 1/(Ngr*(If/Iw-1));                             % Output gain
Kf = 1/Ky;                                          % Setpoint gain
Kb = Ncrev*T_PID*max2absmax/(2*pi);                 % Feedback gain
Kr = Ncrev*T_PID*Ngr*(If/Iw-1)*max2absmax/(2*pi);   % Reference gain
Ku = tau_max/Ndc_max;                               % Input gain
Kh = Ngr*(1-Iw/If)/Iw;                              % Plant gain

K = Kb*Ku*Kh*T_PID;
a1 = -2;
a2 = 1;
b1 = K*Kp;
b2 = K*(Ki*T_PID-Kp);

%% Define LBMPC inputs
N = 30;         % Length of MPC horizon
m = 3;          % Number of total inputs, 3 real and 2 slack 
n = 6;          % Number of states
mslack = 2;     % Number of slack inputs
mtot = m + mslack; % Number of total inputs
T_MPC = 1;      % MPC time step (s)
T_PWM = 0.05;   % PWM step time (s) (Should be divisible into T_MPC)
N_PWM = T_MPC/T_PWM-2;  % Number of PWM steps total
N_PID = T_MPC/T_PID;    % Number of PID steps in an MPC step
DC = N_PWM*T_PWM/T_MPC; % Maximum thruster duty cycle
gamma = 20*pi/180;      % Half angle of LOS cone constraint
rtol = .03;             % Proportion of the radius of the target to subtract
                        % off from the point of the cone. 1 means the cone
                        % is pointed at the center of the target object

% Position state: x_dot y_dot x y
Acxy = [zeros(2) zeros(2)    % State matrix
        eye(2)   zeros(2)];
% Need to add extra states
% Input matrix (two thruster simultaneously firing)
Bcxy = [2*DC*F/mass*eye(2) 
    zeros(2)];             
Ccxy = [zeros(2) eye(2)];    % Output matrix

% Convert to discrete:
sysxy = ss(Acxy,Bcxy,Ccxy,zeros(2));
dsysxy = c2d(sysxy,T_MPC);
Axy = dsysxy.a;
Bxy = dsysxy.b;
Cxy = dsysxy.c;

% Attitude state: theta_dot theta nu (nonphysical parameter)
% Make sure to calibrate axes correctly in Vicon initially, then do a
% coordinate transformation to get relative angle between the two
c = 1;
Ath1 = [-a2*b1/b2-b1  0    % For 1 PID time step
         T_PID        1];
Bth1 = [b1 0]';
Cth1 = [0 1];

Ath = Ath1^N_PID;
Bth = zeros(size(Ath1*Bth1));
for i = 1:N_PID
Bth = Ath1^(N_PID-i)*Bth1+Bth;
end
Bth
Cth = Cth1;

A = [Axy        zeros(4,2)
     zeros(2,4) Ath];
B = [Bxy        zeros(4,3)
     zeros(2,2) Bth    zeros(2,2)];
C = [Cxy        zeros(2,2)
     zeros(1,4) Cth];
D = zeros(3,m+mslack);
k = zeros(n,1);             % Affine offset

Bshort = B(:,1:m);
dsys = ss(A,Bshort,C,D(:,1:m),T_MPC);

% Use LQR to select gains
Q = diag([5000,5000,50,50,50,500]);
Q_tilde = Q;            % p.d. weight matrix for state
R = diag([10,10,1000, 10^5, 10^5]);      % p.s.d. weight matrix on input
dlqr_controlweight = 1; % Scales control weight in LQR
K = -dlqr(A,B,Q,dlqr_controlweight*R);
eig(A+B*K)
P = dare(A+B*K, B, Q, dlqr_controlweight*R);
Q_tilde_f = P;          % p.d. weight matrix for final state
eig(A);

%% Feasibility regions
% Limits in target object frame, limits relative position and velocity
% which isn't exactly what we want, but good enough for now if object is
% not moving very fast. 
x1dot_min = -0.2;       % m/s (true limit is ~0.3) 
x1dot_max = 0.2;        % need at least .15 m/s to activate gripper
x2dot_min = -0.2;
x2dot_max = 0.2;
% x1_min = 0.1;         % meters (true limits)
% x1_max = 3.5;   
% x2_min = 0.1;
% x2_max = 2.5;
x1_min = -3;           % meters (relaxed limits)
x1_max = 3;   
x2_min = -3;
x2_max = 3;
thetadot_min = -pi;     % rad/s (made up)
thetadot_max = pi; 
theta_min = -50*pi;     % rad (in principle, should have no limits)
theta_max = 50*pi;
u1_min = -1;            % Two thrusters, normalized
u1_max = 1;             % Cut actuation limit in half due to 50% actuation duty cycle
u2_min = -1;
u2_max = 1;
u3_min = -20*pi/180;    % [rad/s] Desired rotation rate
u3_max = 20*pi/180;
state_uncert = [1e-3 1e-3 1e-4 1e-4 1e-3 1e-2 1e-3]'; % Made-up numbers

% Set final state constraint to maximally invariant set
model = LTISystem('A',A,'B',Bshort);
model.x.min = [x1dot_min x2dot_min x1_min x2_min thetadot_min theta_min]';
model.x.max = [x1dot_max x2dot_max x1_max x2_max thetadot_max theta_max]';
model.u.min = [u1_min u2_min u3_min];
model.u.max = [u1_max u2_max u3_max];
InvSet = model.invariantSet('maxIterations',40)

h_u = [u1_max u2_max u3_max 10^5 10^5 u1_min u2_min u3_min -10^5 -10^5]';
h_x = [x1dot_max  x2dot_max  x1_max  x2_max  thetadot_max  theta_max...
      x1dot_min x2dot_min x1_min x2_min thetadot_min theta_min]';
h_g = [state_uncert; state_uncert];

