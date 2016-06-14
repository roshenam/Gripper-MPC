function [] = Animate(t, x, XD, rA, rB,rtol, gamma,cone,animate,filename)
if animate
v = VideoWriter(filename);
v.FrameRate = 6;
open(v);
end
x = x';
xrange = linspace(0,10);
phi = XD(3); x0 = x(1,1)+XD(1); y0 = x(1,2)+XD(2); theta0 = x(1,3)+XD(3);
h = figure('Units','Normalized','Position',[.2,.1,.6,.6]);
hold all
xpos = x(:,1)+XD(1); ypos = x(:,2)+XD(2); thetapos = x(:,3)+XD(3);
viscircles([0,0],rA,'color','k','linewidth',1);
dim = max([x(1,1)+XD(1),x(1,2)+XD(2)]);
xlim([-2,dim+1])
ylim([-2,dim+1])
axis('square')
viscircles([rA*cos(phi), rA*sin(phi)],rA/10,'color','k');
y1 = (sin(phi+gamma).*xrange./((rA-rtol)*sin(gamma))-1)./(cos(phi+gamma)/((rA-rtol)*sin(gamma)));
y2 = (1+ sin(phi-gamma).*xrange./((rA-rtol)*sin(gamma)))./(cos(phi-gamma)/((rA-rtol)*sin(gamma)));
if cone
plot(xrange,y1,'k','linewidth',2)
plot(xrange,y2,'k','linewidth',2)
end
freeflyer = rectangle('Position',[x0-rB, y0-rB, rB*2, rB*2],'Curvature',[1,1],'facecolor','g');
gripper = rectangle('Position',[x0-rB*cos(theta0)-rB/7, y0-rB*sin(theta0)-rB/7,2*rB/7,2*rB/7],'Curvature',[1,1],'facecolor','k');
path = line(x(1,1),x(1,2),'color','b');

frame = getframe;
if animate
writeVideo(v,frame);
end
for i=2:length(x)
    set(freeflyer, 'Position', [xpos(i)-rB, ypos(i)-rB, rB*2, rB*2]);
    set(gripper, 'Position', [xpos(i)-rB*cos(thetapos(i))-rB/7, ypos(i)-rB*sin(thetapos(i))-rB/7,2*rB/7,2*rB/7]);
    set(path,'XData',x(1:i,1),'YData',x(1:i,2))
    drawnow  
    pause(.1)
    if animate
    frame = getframe;
    writeVideo(v,frame);
    end
end
if animate
close(v);
end