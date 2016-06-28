function [] = Animate(init, params, xtot, animate, filename, rotate)
% If animate is 1, a VideoWriter object will be created and a video will be
% saved of the animation. If rotate is 1, it will plot the rotating cone
% constraint. xtot must include phi data. 
if animate
v = VideoWriter(filename);
v.FrameRate = 6;
open(v);
end
xtot = xtot';
xrange = linspace(-4,10);
rA = params.rp; rtol = params.rtol; rB = params.rs;
gamma = params.gamma;
phi = params.phi; x0 = xtot(1,1); y0 = xtot(1,2); theta0 = xtot(1,3);
h = figure('Units','Normalized','Position',[.2,.1,.6,.6]);
hold all
xpos = xtot(:,1); ypos = xtot(:,2); thetapos = xtot(:,3)+phi;
viscircles([0,0],rA,'color','k','linewidth',1);
dim = max([max(xtot(:,1)),max(xtot(:,2))]);
xlim([-4,dim+1])
ylim([-4,dim+1])
axis('square')
viscircles([rA*cos(phi), rA*sin(phi)],rA/10,'color','k');
y1 = (sin(phi+gamma).*xrange./((rA-rtol)*sin(gamma))-1)./(cos(phi+gamma)/((rA-rtol)*sin(gamma)));
y2 = (1+ sin(phi-gamma).*xrange./((rA-rtol)*sin(gamma)))./(cos(phi-gamma)/((rA-rtol)*sin(gamma)));
cone1 = line(xrange,y1,'color','k','linewidth',2);
cone2 = line(xrange,y2,'color','k','linewidth',2);

freeflyer = rectangle('Position',[x0-rB, y0-rB, rB*2, rB*2],'Curvature',[1,1],'facecolor','g');
%gripper = rectangle('Position',[x0-rB*cos(theta0)-rB/7, y0-rB*sin(theta0)-rB/7,2*rB/7,2*rB/7],'Curvature',[1,1],'facecolor','k');
path = line(xtot(1,1),xtot(1,2),'color','b');

frame = getframe;
if animate
writeVideo(v,frame);
end
for i=2:length(xtot)
    set(freeflyer, 'Position', [xpos(i)-rB, ypos(i)-rB, rB*2, rB*2]);
    %set(gripper, 'Position', [xpos(i)-rB*cos(thetapos(i))-rB/7, ypos(i)-rB*sin(thetapos(i))-rB/7,2*rB/7,2*rB/7]);
    set(path,'XData',xtot(1:i,1),'YData',xtot(1:i,2))
    if rotate
        set(cone1,'YData',(sin(xtot(i,end)+gamma).*xrange./((rA-rtol)*sin(gamma))-1)./(cos(xtot(i,end)+gamma)/((rA-rtol)*sin(gamma))));
        set(cone2,'Ydata',(1+ sin(xtot(i,end)-gamma).*xrange./((rA-rtol)*sin(gamma)))./(cos(xtot(i,end)-gamma)/((rA-rtol)*sin(gamma))));
    end
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