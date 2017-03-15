% Defining parameters

clear all

phi = pi/4; omega = 2*pi/180;
params.phi = phi; % Phi is the angle of the port from the defined horizontal
params.omega = omega; % Omega is the rotation in rad/sec of the target
params.rp = .15; params.rtol = params.rp+.02;  params.rs = .15; params.gamma = 50*pi/180;
params.w = .01;
params.l = .1;
%theta_c = 4*pi/180;
theta_c = 0*pi/180;
gamma_max = atan(.05/(params.rp+params.rs)) - theta_c
%gamma_max = 5*pi/180;
params.gamma = gamma_max;
xs = 6;
ys = 1;
%for i = 1:100
gamma0 = atan(.05/(params.rp+params.rs));
omega_c = 100*pi/180;
% rp is the radius of the target, rs is the radius of the spacecraft, gamma
% is the half angle of the constraint cone
params.Umax = .2;
% Max thrust in Newtons/kg of the free flyer
params.Tmax = .3; % 50% of stall torque
% Max torque in N*m/(kg*m^2)
params.Ts = 0.5;
% Discretization constant (in seconds)
params.Qval = 10^3; params.Rval = 10^2;
params.Qval2 = 10^2;
% Weights on states, control inputs, and slack variables
params.betaHIGH = -1.5; params.betaLOW = -0.2;
params.eta = 1;
params.tantol = .4;
% Constants defining the shape of the bounding function for normal
% component of hte velocity

rp = params.rp; rs = params.rs;
Ts = params.Ts;
betaHIGH = params.betaHIGH; betaLOW = params.betaLOW;
eta = params.eta; tantol = params.tantol;
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

%Dd = zeros(5,3);
%Dd = zeros(3,3);
Dd = zeros(12,5);
Cd = zeros(12,6);
%Cd = zeros(2,6);
eta_vel = vmax/(sqrt(2)*(rp+rs));
eta_omega = omega_c/(rp+rs);

%a = 0;
%a = 10;
%b = gamma0;
%b = 30;
gamma = gamma0;
N_vec = zeros(length(ys), length(xs));
%[X,Y] = meshgrid(x,y);
counter = 0;
for i = 1:length(xs)
    for j = 1:length(ys)
        %counter = counter + 1;
        x = xs(i)
        y = ys(j)
        %b = gamma0;
        %a = 0;
        b = 30;
        a = 10;
        params.N = b;
        %a = 10;
        %gamma = b;
        %params.gamma = b;
        error_flag = 0;
        cprev = nan;
        c = nan;
        check = 1;
        % Check if feasible at limit b
        check_feasible
        if check == -1          % If not feasible at limit b, stop
            gamma_vec = nan;
        else                    % If feasible at limit b, do bisection
            while b-a >= 5;
                error_flag = 0;
                disp(['a = ', num2str(a)])
                disp(['b = ', num2str(b)])
                %cprev = c;
                c = floor(a+(b-a)/2)
                %params.gamma = c;
                params.N = c;
                %gamma = c;
                check = 1;
                check_feasible
                if check == -1
                    a = c;
                    cprev = b;
                else
                    b = c;
                    cprev = c;
                end

                
            end
            %check
        end
        %gamma_vec(counter) = c;
        %N_vec(counter) = cprev;
        N_vec(j, i) = cprev;
    end
end




