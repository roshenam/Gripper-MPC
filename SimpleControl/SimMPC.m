function [times, states] = SymMPC(ParamsS,ParamsT,Ts1,Ts2,N,tol,Tset,animate)
% Roshena MacPherson: March 30, 2017
%
% ParamsS must contain the radius of the spacecraft rs, umax its maximum thrust in N/kg
% tmax its maximum torque in N*m/(kg*m^2), xinit which is a 6 element column vector
% containing the spacecraft's initial conditions specified in the table
% frame, vf which is the final velocity we want to hit the target with
%
% ParamsT must contain the radius of the target rt, omega its angular velocity in
% the table frame, nu0 the initial offset angle of the port from the table
% horizontal
%
% Ts1 is the MPC discretization time in phase 1
%
% Ts2 is the discretization time in phase 2 (can be shorter because no
% optimization is happening)
%
% N is the horizon length for the MPC
%
% tol is the tolerance for how aligned the spacecraft's velocity vector is
% with its position vector before we start the second phase
%
% Tset is the terminal set for the initial MPC phase computed using
% CalcTSet
%
% if animate is set to 1 the result is animated
%% Setting up system
A = zeros(4,4);
A(1,3) = 1; A(2,4) = 1;
B = zeros(4,2);
B(3,1) = 1; B(4,2) = 1;
C = zeros(2,4);
C(1,1) = 1; C(2,2) = 1;
D = zeros(2,2);
sysC = ss(A,B,C,D);
sysD = c2d(sysC,Ts1);
Ad = sysD.A; Bd = sysD.B;
system = LTISystem('A',Ad,'B',Bd,'Ts',Ts1);
umax = ParamsS.umax;
Q = eye(4); Q(3,3) = 100; Q(4,4) = 100;
R = 10^2.*eye(2);
system.x.min = [-2; -2; -1; -1]; %to stay on the table, max speed allowed on table
system.x.max = [2; 2; 1; 1];
system.u.min = [-umax; -umax];
system.u.max = [umax; umax];
system.x.penalty = QuadFunction(Q);
system.u.penalty = QuadFunction(R);
Pn = system.LQRPenalty;
%Tset = system.LQRSet;
%Tset = system.invariantSet();
system.x.with('terminalPenalty');
system.x.terminalPenalty = Pn;
system.x.with('terminalSet');
system.x.terminalSet = Tset;
mpc = MPCController(system,N);
loop = ClosedLoop(mpc, system);

%% Simulating alignment portion
Nsim = 60;
x0vec = ParamsS.xinit(1:4);
data = loop.simulate(x0vec, Nsim)
xout = data.X; xout = xout';
rvec = zeros(Nsim+1,3); vvec = zeros(Nsim+1,3);
rvec(:,1:2) = xout(:,1:2); vvec(:,1:2) = xout(:,3:4);
error = cross(rvec,vvec); %we want the velocity and position vectors to be parallel for alignment
error = error(:,3);
%plot(1:(Nsim+1),error,'ro')
%plot(xout(:,1),xout(:,2))
idx = find(abs(error)<tol,1)

%% Setting up final approach
% for the second phase we're going to rotate into the frame with the
% spacecraft on the negative x axis
phi = atan2(xout(idx,2),xout(idx,1)) + pi  
Rn2r = [cos(phi) sin(phi); -sin(phi) cos(phi)];
Rr2n = Rn2r';
xout_rot = zeros(idx,4);
% rotate position data into the new frame
for i=1:idx
    xout_rot(i,1:2) = Rn2r*xout(i,1:2)';
    xout_rot(i,3:4) = Rn2r*xout(i,3:4)';
end

% %%
% rt = .15; rs = .15;
% h = figure('Units','Normalized','Position',[.2,.1,.6,.6]);
% spacecraft = rectangle('Position',[xout(1,1)-rs, xout(1,2)-rs, rs*2, rs*2],'Curvature',[1,1],'facecolor',[0 1 1]);
% target = rectangle('Position',[-rt, -rt, 2*rt, 2*rt],'Curvature',[1,1],'facecolor',[1 1 0]);
% 
% xlim([-4,4])
% xsize = 8;
% ylim([-xsize/2, xsize/2])
% axis('square')
% path = line(xout(1,1),xout(1,2),'color','b','linewidth',2);
% 
% for i=2:idx
%     set(spacecraft, 'Position', [xout(i,1)-rs, xout(i,2)-rs, rs*2, rs*2]);
%     set(path,'XData',xout(1:i,1),'YData',xout(1:i,2))
%     drawnow
%     pause(.1)
% end
%% Calculate control effort for final approach
% Assume here we're working in frame centered at the target COM
xt = 0; rt = ParamsT.rt; rs = ParamsS.rs;
omega = ParamsT.omega;
nu0 = ParamsT.nu0;
% Need the offset angle of the target location from the negative horizontal
nu_n = -(omega/abs(omega))*wrapTo2Pi((phi-pi-wrapTo2Pi(nu0 + omega*Ts1*(idx))))
vf = ParamsS.vf;
x0 = xout_rot(end,1) 
% for this, assuming our velocity is only along the horizontal
% not true in general, but will use LQR to force us towards this trajectory
v0 = xout_rot(end,3) 
xf = (xt-rt-rs); % final location of center of spacecraft
Rdis = xf-x0
t1 = abs((vf-v0)/umax) % time to accelerate/deccelerate from v0 to vf at max thrust
if(vf>v0)
    d1 = .5*umax*t1^2+v0*t1; % distance traveled over t1
else
    d1 = .5*(-umax)*t1^2+v0*t1;
end
k = (abs(omega)*t1-nu_n+abs(omega)*(Rdis-d1)/vf)/(2*pi);
if(vf>v0) % depending on if accel or deccel round up or down
    N = ceil(k);
else
    N = floor(k);
end

twait = (N-k)*2*pi/(abs(omega)*(1-v0/vf)) % initial coasting wait time for alignment
t2 = (Rdis-d1-twait*v0)/vf % coasting time
if(t2<0)
    disp('Not possible!! Move away and try again')
    xout = NaN;
    tout = NaN;
    return;
end
n_phase2 = floor((t1+t2+twait)/Ts2) %number of simulation steps
t = linspace(0,t1+t2+twait,n_phase2);
u_desrot = zeros(length(t),2); % vector of desired thrust values in rotated frame
idx1 = find(t>twait,1);
idx2 = find(t>(twait+t1),1);

if(vf<v0)
    u_desrot(idx1:idx2,1) = -umax;
else
    u_desrot(idx1:idx2,1) = umax;
end
%plot(u_desrot(:,1),'ro')
u_des = zeros(size(u_desrot));
% rotate back to initial frame
for i=1:n_phase2
    u_des(i,:) = Rr2n*u_desrot(i,:)'; % rotating thrust back to table frame
end

%% Generate states for desired trajectory
x_desrot = zeros(length(u_des),2);
v_desrot = zeros(length(u_des),2);
taccel = t(idx1:idx2);
tcoast = t(idx2:end);
x_desrot(1:idx1,1) = x0+v0.*t(1:idx1);
v_desrot(1:idx1,1) = v0;
if(vf<v0)
    x_desrot(idx1:idx2,1) = x_desrot(idx1,1)+v0.*(taccel-twait)-.5*umax.*((taccel-twait).^2);
    v_desrot(idx1:idx2,1) = v0-umax.*(taccel-twait);
else
    x_desrot(idx1:idx2,1) = x_desrot(idx1,1)+v0.*(taccel-twait)+.5*umax.*((taccel-twait).^2);
    v_desrot(idx1:idx2,1) = v0+umax.*(taccel-twait);
end
x_desrot(idx2:end,1) = x_desrot(idx2,1) + vf.*(tcoast-twait-t1);
v_desrot(idx2:end,1) = vf;
% rotate back to table frame
x_des = zeros(size(x_desrot));
v_des = zeros(size(v_desrot));
for i=1:length(x_des)
    x_des(i,:) = Rr2n*x_desrot(i,:)';
    v_des(i,:) = Rr2n*v_desrot(i,:)';
end

state_des = [x_des v_des]; % vector of desired states for phase 2 trajectory
% figure
% plot(x_des(:,1),x_des(:,2))
%% Simulate phase 2
[K,~,~] = lqr(sysC,Q,R); %Compute LQR controller gains
y0 = xout(idx,:);
tspan = [idx*Ts1,idx*Ts1+t1+t2+twait];
Pstruct.A = A;
Pstruct.B = B;
Pstruct.u = u_des;
Pstruct.t = t+idx*Ts1;
Pstruct.i = 1;
Pstruct.K = K;
Pstruct.umax = umax;
Pstruct.state_des = state_des;
%figure(20)
[tout,yout] = ode45(@(t,x) spacecraft_dyn(t,x,Pstruct),tspan,y0);
%plot(yout(:,1),yout(:,2))
%ylim([0,2])
%xlim([-1,1])
%axis('square')

%xtotal = [xout(1:idx,1); yout(1:end,1)];
%ytotal = [xout(1:idx,2); yout(1:end,2)];
xtotal = [xout(1:idx,1); yout(2:end,1)];
ytotal = [xout(1:idx,2); yout(2:end,2)];

%%

%thetas = theta + omega.*(tout-tout(1))-phi;
%theta_first = thetas(1) - omega*Nsim*idx;
t_phase1 = linspace(0,idx*Ts1,idx);
t_phase2 = tout(2:end)';
%theta_firstvec = theta_first + omega.*t_phase1;
nu_phase1 = nu0 + omega.*t_phase1;
nu_phase2 = nu_phase1(end) + omega.*(t_phase2-t_phase2(1)+(t_phase2(2)-t_phase2(1)));
nu_total = [nu_phase1 nu_phase2]';
%thetatotal = [ theta_firstvec thetas' ];
%dockpos = [-rt.*cos(thetatotal); rt.*sin(thetatotal)];
%dockpos_rot = [-rt.*cos(thetas); rt.*sin(thetas)];
dockpos = [-rt.*cos(nu_total) rt.*sin(nu_total)];

ttotal = zeros(length(xtotal),1);
ttotal(1:idx) = t_phase1;
ttotal((idx+1):end) = t_phase2;
%%
if animate==1
    %h = figure('Units','Normalized','Position',[.2,.1,.6,.6]);
    figure(3)
    spacecraft = rectangle('Position',[xtotal(1)-rs, ytotal(1)-rs, rs*2, rs*2],'Curvature',[1,1],'facecolor',[0 1 1]);
    target = rectangle('Position',[-rt, -rt, rt*2, rt*2],'Curvature',[1,1],'facecolor',[1 1 0]);
    xlim([-4,4])
    xsize = 8;
    ylim([-xsize/2, xsize/2])
    axis('square')
    port = rectangle('Position',[dockpos(1,1)-rt/8,dockpos(1,2)-rt/8, rt/4, rt/4],'Curvature',[1,1],'facecolor','r');
    path = line(xtotal(1),ytotal(1),'color','g','linewidth',1);
    for i=2:length(dockpos)
        set(spacecraft, 'Position', [xtotal(i)-rt, ytotal(i)-rt, rt*2, rt*2]);
        set(port, 'Position', [dockpos(i,1)-rt/8,dockpos(i,2)-rt/8, rt/4, rt/4]);
        set(path,'XData',xtotal(1:i),'YData',ytotal(1:i))
        drawnow
        %pause((ttotal(i)-ttotal(i-1))/10)
    end
elseif animate==2
        %h = figure('Units','Normalized','Position',[.2,.1,.6,.6]);
    figure(3)
    spacecraft = rectangle('Position',[xtotal(end)-rs, ytotal(end)-rs, rs*2, rs*2],'Curvature',[1,1],'facecolor',[0 1 1]);
    target = rectangle('Position',[-rt, -rt, rt*2, rt*2],'Curvature',[1,1],'facecolor',[1 1 0]);
    xlim([-4,4])
    xsize = 8;
    ylim([-xsize/2, xsize/2])
    axis('square')
    port = rectangle('Position',[dockpos(end,1)-rt/8,dockpos(end,2)-rt/8, rt/4, rt/4],'Curvature',[1,1],'facecolor','r');
    path = line(xtotal(1),ytotal(1),'color','g','linewidth',1);
    set(path,'XData',xtotal(1:end),'YData',ytotal(1:end))
    drawnow
end

times = ttotal;
states = [xtotal ytotal nu_total];


