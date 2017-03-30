A = [0 1; 0 0]; B = [0 1]';
C = eye(2); D = zeros(2,1);
Ts = .5;
sysC = ss(A,B,C,D);
sysD = c2d(sysC, Ts);
Ad = sysD.A;
Bd = sysD.B;

xt = 0; rt = .15; rs = .15;
omega = -400*pi/180;
theta = pi/4;
x0 = -0.8586; v0 = 0.0854; xf = (xt-rt-rs); vf = .2;
R = xf-x0
umax = .4/18;
t1 = abs((vf-v0)/umax);
if(vf>v0)
    d1 = .5*umax*t1^2+v0*t1
else
    d1 = .5*(-umax)*t1^2+v0*t1
end
%d2 = R-d1;
%t2 = d2/vf;
k = (abs(omega)*t1-theta+abs(omega)*(R-d1)/vf)/(2*pi)
if(vf>v0)
    N = ceil(k)
else
    N = floor(k)
end

twait = (N-k)*2*pi/(abs(omega)*(1-v0/vf))
t2 = (R-d1-twait*v0)/vf
t1
t = linspace(0,t1+t2+twait,500);
xs = zeros(length(t),1);
idx1 = find(t>twait,1);
idx2 = find(t>(twait+t1),1);
taccel = t(idx1:idx2);
tcoast = t(idx2:end);
xs(1:idx1) = x0+v0.*t(1:idx1);
if(vf<v0)
    xs(idx1:idx2) = xs(idx1)+v0.*(taccel-twait)-.5*umax.*((taccel-twait).^2);
else
    xs(idx1:idx2) = xs(idx1)+v0.*(taccel-twait)+.5*umax.*((taccel-twait).^2);
end
xs(idx2:end) = xs(idx2) + vf.*(tcoast-twait-t1);
plot(t,xs)

%%
thetas = theta + omega.*t;
dockposx = xt-rt.*cos(thetas);
dockposy = rt.*sin(thetas);
h = figure('Units','Normalized','Position',[.2,.1,.6,.6]);
spacecraft = rectangle('Position',[x0-rs, -rs, rs*2, rs*2],'Curvature',[1,1],'facecolor',[0 1 1]);
target = rectangle('Position',[xt-rt, -rt, rt*2, rt*2],'Curvature',[1,1],'facecolor',[1 1 0]);
xlim([x0-2,xt+3])
xsize = xt+2-(x0-2);
ylim([-xsize/2, xsize/2])
axis('square')
port = rectangle('Position',[dockposx(1)-rt/8,dockposy(1)-rt/8, rt/4, rt/4],'Curvature',[1,1],'facecolor','r');

for i=2:length(xs)
    set(spacecraft, 'Position', [xs(i)-rt, -rt, rt*2, rt*2]);
    set(port, 'Position', [dockposx(i)-rt/8,dockposy(i)-rt/8, rt/4, rt/4]);
    drawnow  
    pause(t(2))
end