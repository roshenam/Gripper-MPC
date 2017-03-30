%% Setting up system
A = zeros(4,4);
A(1,3) = 1; A(2,4) = 1;
B = zeros(4,2);
B(3,1) = 1; B(4,2) = 1;
C = zeros(2,4);
C(1,1) = 1; C(2,2) = 1;
D = zeros(2,2);
sysC = ss(A,B,C,D);
Ts = .5;
sysD = c2d(sysC,Ts);
Ad = sysD.A; Bd = sysD.B; Cd = sysD.C; Dd = sysD.D;
system = LTISystem('A',Ad,'B',Bd,'Ts',Ts);
umax = .4/18;
Q = eye(4); Q(3,3) = 100; Q(4,4) = 100;
R = 10^2.*eye(2);
system.x.min = [-4; -4; -1; -1];
system.x.max = [4; 4; 1; 1];
system.u.min = [-umax; -umax];
system.u.max = [umax; umax];
system.x.penalty = QuadFunction(Q);
system.u.penalty = QuadFunction(R);
Pn = system.LQRPenalty;
%Tset = system.LQRSet;
Tset = system.invariantSet();
system.x.with('terminalPenalty');
system.x.terminalPenalty = Pn;
system.x.with('terminalSet');
system.x.terminalSet = Tset;
N = 20;
mpc = MPCController(system,N);
loop = ClosedLoop(mpc, system);

%% Simulating alignment portion
Nsim = 60;
x0vec = [1,1,.2,.3]';
data = loop.simulate(x0vec, Nsim)
xout = data.X; xout = xout';
atan2(xout(1,2),xout(1,1))-atan2(xout(1,4),xout(1,3))
%error = atan(xout(:,2)./xout(:,1))-atan(xout(:,4)./xout(:,3));
rvec = zeros(Nsim+1,3); vvec = zeros(Nsim+1,3);
rvec(:,1:2) = xout(:,1:2); vvec(:,1:2) = xout(:,3:4);
error = cross(rvec,vvec);
error = error(:,3);
plot(1:(Nsim+1),error,'ro')

idx = find(abs(error)<.001,1)

%% Setting up final approach
phi = atan2(xout(idx,2),xout(idx,1)) + pi
Rn2r = [cos(phi) sin(phi); -sin(phi) cos(phi)];
Rr2n = Rn2r';
xout_rot = zeros(idx,4);
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
Ts_phase2 = .01;
xt = 0; rt = .15; rs = .15;
omega = -400*pi/180;
theta = pi/4; % change this eventually
x0 = xout_rot(end,1) 
v0 = xout_rot(end,3) 
xf = (xt-rt-rs); vf = .2;
R = xf-x0
umax = .4/18;
t1 = abs((vf-v0)/umax)
if(vf>v0)
    d1 = .5*umax*t1^2+v0*t1;
else
    d1 = .5*(-umax)*t1^2+v0*t1;
end
%d2 = R-d1;
%t2 = d2/vf;
k = (abs(omega)*t1-theta+abs(omega)*(R-d1)/vf)/(2*pi);
if(vf>v0)
    N = ceil(k);
else
    N = floor(k);
end

twait = (N-k)*2*pi/(abs(omega)*(1-v0/vf))
t2 = (R-d1-twait*v0)/vf
n_phase2 = floor((t1+t2+twait)/Ts_phase2) %number of simulation steps
t = linspace(0,t1+t2+twait,n_phase2);
u_desrot = zeros(length(t),2);
idx1 = find(t>twait,1);
idx2 = find(t>(twait+t1),1);

if(vf<v0)
    u_desrot(idx1:idx2,1) = -umax;
else
    u_desrot(idx1:idx2,1) = umax;
end
plot(u_desrot(:,1),'ro')
u_des = zeros(size(u_desrot));
% rotate back to initial frame
for i=1:n_phase2
    u_des(i,:) = Rr2n*u_desrot(i,:)';
end

%% Simulate phase 2
y0 = xout(idx,:);
tspan = [idx*Ts,idx*Ts+t1+t2+twait];
Pstruct.A = A;
Pstruct.B = B;
Pstruct.u = u_des;
Pstruct.t = t+idx*Ts;
Pstruct.i = 1;
[tout,yout] = ode45(@(t,x) spacecraft_dyn(t,x,Pstruct),tspan,y0);
plot(yout(:,1),yout(:,2))
ylim([0,2])
xlim([-1,1])
axis('square')
%%
xtotal = [xout(1:idx,1); yout(1:end,1)];
ytotal = [xout(1:idx,2); yout(1:end,2)];


%%

thetas = theta + omega.*(tout-tout(1))-phi;
theta_first = thetas(1) - omega*Nsim*idx;
t_phase1 = linspace(0,idx*Ts,idx);
theta_firstvec = theta_first + omega.*t_phase1;
thetatotal = [ theta_firstvec thetas' ];
dockpos = [-rt.*cos(thetatotal); rt.*sin(thetatotal)];
dockpos_rot = [-rt.*cos(thetas); rt.*sin(thetas)];

%%
h = figure('Units','Normalized','Position',[.2,.1,.6,.6]);
spacecraft = rectangle('Position',[xtotal(1)-rs, ytotal(1)-rs, rs*2, rs*2],'Curvature',[1,1],'facecolor',[0 1 1]);
target = rectangle('Position',[-rt, -rt, rt*2, rt*2],'Curvature',[1,1],'facecolor',[1 1 0]);
xlim([-4,4])
xsize = 8;
ylim([-xsize/2, xsize/2])
axis('square')
port = rectangle('Position',[dockpos(1,1)-rt/8,dockpos(2,1)-rt/8, rt/4, rt/4],'Curvature',[1,1],'facecolor','r');
path = line(xtotal(1),ytotal(1),'color','b','linewidth',1);
ttotal = zeros(length(xtotal),1);
ttotal(1:idx) = linspace(0,idx*Ts,idx);
ttotal((idx+1):end) = tout;
for i=2:length(dockpos)
    set(spacecraft, 'Position', [xtotal(i)-rt, ytotal(i)-rt, rt*2, rt*2]);
    set(port, 'Position', [dockpos(1,i)-rt/8,dockpos(2,i)-rt/8, rt/4, rt/4]);
    set(path,'XData',xtotal(1:i),'YData',ytotal(1:i))
    drawnow 
    pause((ttotal(i)-ttotal(i-1))/10)  
end




