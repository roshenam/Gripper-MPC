% Gen_InvSet
% This function takes as input a structural array of parameters including
% radii, mass, thrust, as well as the angular velocity of the target object
% in rad/sec. It outputs a matrix A and vector b defining the polyhedral
% terminal set for the system. errorflag is 1 if the sets are empty
function [A, b, errorflag] = Gen_InvSet( params, omega )

% Constants defining the shape of the bounding function for normal
% component of hte velocity
gamma = atan(params.dmax/(params.rs+params.w)) - params.theta_c;
rt = params.rt; rs = params.rs;
Ts = params.Ts;
vmax = .3; % velocity used on d and omega plots from old paper
% Creating dynamic system
A = zeros(6,6); A(1,4) = 1; A(2,5) = 1; A(3,6) = 1;
A(4,1) = omega^2; A(4,5) = 2*omega;
A(5,2) = omega^2; A(5,4) = -2*omega;
B = zeros(6,5); B(4,1) = 1; B(5,2) = 1; B(6,3) = 1;
C = zeros(6,6);  D = zeros(6,5);
sysD = ss(A,B,C,D);
sysD = c2d(sysD,Ts);

Ad = sysD.a; Bd = sysD.b; % Matrices for discrete system

% Y is a vector of constrained outputs. C and D are matrices defining Y's
% dependence on the states and control inputs

Dd = zeros(10,5);
Cd = zeros(10,6);
eta_vel = vmax/(sqrt(2)*(rt+rs));
Cd(1,1) = -tan(pi/4+gamma); Cd(1,2) = 1; Dd(1,4) = 1; % LOS constraints
Cd(2,1) = tan(pi/4-gamma); Cd(2,2) = -1; Dd(2,5) = 1;
Cd(3,1) = -eta_vel; Cd(3,2) = -eta_vel; Cd(3,4) = -1; % Velocity constraints
Cd(4,1) = -eta_vel; Cd(4,2) = -eta_vel; Cd(4,4) = 1;
Cd(5,1) = -eta_vel; Cd(5,2) = -eta_vel; Cd(5,5) = 1; 
Cd(6,1) = -eta_vel; Cd(6,2) = -eta_vel; Cd(6,5) = -1;
Cd(7,4) = 1; % Angle of attack constraints through velocity 
Cd(8,5) = 1; 
Cd(9,1) = -1; % x > 0
Cd(10,2) = -1; % y > 0 

% Upper bounds on thrust inputs and slack variables
Umax = params.Umax;
Tmax = params.Tmax; % max torque applied by flywheel
Ymax = [(params.rt-params.rtol)/(1-tan(pi/4+gamma)); -(params.rt-params.rtol)/(1-tan(pi/4-gamma)); 0; 0; 0; 0; 0; 0; 0; 0];%; 100*vmax/sqrt(2)];...

% Creating weighting matrices for cost function
Q = params.Qval.*eye(6); Q(3,3) = 10000;
R = params.Rval.*eye(3); R(3,3) = 10; R(4,4) = params.slackweight; R(5,5) = params.slackweight;
system = LTISystem('A',Ad,'B',Bd,'C',Cd,'D',Dd,'Ts',Ts);

system.u.min = [-Umax; -Umax; -Tmax; -10^4; -10^4];
system.u.max = [Umax; Umax; Tmax; 10^4; 10^4];
system.y.max = Ymax;
system.x.penalty = QuadFunction(Q);
system.u.penalty = QuadFunction(R);
Pn = system.LQRPenalty;
system.x.with('terminalPenalty');
system.x.terminalPenalty = Pn;
Tset = system.LQRSet;
A = Tset.A;
b = Tset.b;

if( isempty(A))
    errorflag = 1;
    disp('The set is empty');
else
    errorflag = 0;
    
end