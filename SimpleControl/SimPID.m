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
Ad = sysD.A;
Bd = sysD.B;
Ax = [0 1; 0 0]; Bx = [0 1]';
P = [-1; -2];
K = acker(Ax,Bx,P)

%%
x0vec = [1,1.5,-.1,0]';
kp = .001;
ki = .001;
kd = 0;
xout = zeros(1,4);
xout(1,:) = x0vec;
uout = zeros(1,2);
xcurr = x0vec;
xlast = x0vec;
integralx = 0;
integraly = 0;
umax = .4/18;
error = abs(x0vec(2)*x0vec(3) - x0vec(1)*x0vec(4))
errorx = x0vec(1);
errory = x0vec(2);
i = 1;
while (error>.01)
%for i=1:400
    x = xcurr(1); y = xcurr(2); vx = xcurr(3); vy = xcurr(4);
    ux = K*[x; vx];
    uy = K*[y; vy];
    if(ux>umax)
        ux = umax;
    elseif(ux < -umax)
        ux = -umax;
    end
    if(uy>umax)
        uy = umax;
    elseif(uy < -umax)
        uy = -umax;
    end
    uvec = [ux; uy];
    xlast = xcurr;
    xcurr = Ad*xcurr - Bd*uvec;
    xout(i+1,:) = xcurr;
    uout(i,:) = uvec;
    i=i+1;
    error = abs(xcurr(2)*xcurr(3) - xcurr(1)*xcurr(4));
end
n=i-1;
%n=400;
subplot(4,1,1)
plot(xout(:,1),xout(:,2))
title('path')
subplot(4,1,2)
plot(1:n,uout(:,1))
title('control effort')
ylim([-umax-.01,umax+.01])
subplot(4,1,3)
plot(1:(n+1),(xout(:,2).*xout(:,3) - xout(:,1).*xout(:,4)))
title('error')
subplot(4,1,4)
%plot(1:(n+1), atan(xout(:,2)./xout(:,1)))
title('atans')
%plot(1:(n+1), atan(xout(:,4)./xout(:,3)))
plot(1:(n+1),xout(:,1))
hold all
plot(1:(n+1),xout(:,2))

rt = .1; rs = .1;
h = figure('Units','Normalized','Position',[.2,.1,.6,.6]);
spacecraft = rectangle('Position',[x0vec(1)-rs, x0vec(2)-rs, rs*2, rs*2],'Curvature',[1,1],'facecolor',[0 1 1]);
target = rectangle('Position',[-rt, -rt, rt*2, rt*2],'Curvature',[1,1],'facecolor',[1 1 0]);

xlim([-2,2])
xsize = 4;
ylim([-xsize/2, xsize/2])
axis('square')
path = line(xout(1,1),xout(1,2),'color','b','linewidth',2);

for i=2:length(xout)
    set(spacecraft, 'Position', [xout(i,1)-rs, xout(i,2)-rs, rs*2, rs*2]);
    set(path,'XData',xout(1:i,1),'YData',xout(1:i,2))
    drawnow  
    pause(.1)
end