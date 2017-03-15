function [] = Animate(init, params, xtot, phidata, animate, filename, rotate, attitude)
% If animate is 1, a VideoWriter object will be created and a video will be
% saved of the animation. If rotate is 1, it will plot the rotating cone
% constraint. xtot must include phi data. If attitude is 1 then xtot should
% have orientation data in it. 
if animate
v = VideoWriter(filename);
v.FrameRate = 6;
open(v);
end
xtot = xtot';
rA = params.rp; rtol = params.rtol; rB = params.rs;
gamma = params.gamma;
phi = params.phi; x0 = xtot(1,1); y0 = xtot(1,2); 
if attitude
theta0 = xtot(1,3);
end
h = figure('Units','Normalized','Position',[.2,.1,.6,.6]);
hold all
xpos = xtot(:,1); ypos = xtot(:,2); 
if attitude
thetapos = xtot(:,3);
end
viscircles([0,0],rA,'color','k','linewidth',1);
dim = max([max(xtot(:,1)),max(xtot(:,2))]);
xrange = linspace(0,dim+1);
xlim([-2,dim+1])
ylim([-2,dim+1])
%xlim([-2,dim])
%ylim([-1,dim+1])
axis('square')
set(gca, 'Ticklength', [0 0])
set(gca, 'XTick',[])
set(gca, 'YTick',[])

port = rectangle('Position',[rA*cos(phi)-rA/8,rA*sin(phi)-rA/8, rA/4, rA/4],'Curvature',[1,1],'facecolor','r');
y1 = (sin(phi+gamma).*xrange./((rA-rtol)*sin(gamma))-1)./(cos(phi+gamma)/((rA-rtol)*sin(gamma)));
y2 = (1+sin(phi-gamma).*xrange./((rA-rtol)*sin(gamma)))./(cos(phi-gamma)/((rA-rtol)*sin(gamma)));
cone1 = line(xrange,y1,'color','k','linewidth',2);
cone2 = line(xrange,y2,'color','k','linewidth',2);

freeflyer = rectangle('Position',[x0-rB, y0-rB, rB*2, rB*2],'Curvature',[1,1],'facecolor',[0 1 1]);
if attitude
gripper = rectangle('Position',[x0-rB*cos(theta0)-rB/4, y0-rB*sin(theta0)-rB/4,2*rB/4,2*rB/4],'Curvature',[1,1],'facecolor','r');
else
end
path = line(xtot(1,1),xtot(1,2),'color','b','linewidth',2);

frame = getframe;
if animate
writeVideo(v,frame);
end
for i=2:length(xtot)
    set(freeflyer, 'Position', [xpos(i)-rB, ypos(i)-rB, rB*2, rB*2]);
    if attitude
    set(gripper, 'Position', [xpos(i)-rB*cos(thetapos(i))-rB/4, ypos(i)-rB*sin(thetapos(i))-rB/4,2*rB/4,2*rB/4]);
    else
    end
    set(path,'XData',xtot(1:i,1),'YData',xtot(1:i,2))
    if rotate
        set(cone1,'YData',(sin(phidata(i)+gamma).*xrange./((rA-rtol)*sin(gamma))-1)./(cos(phidata(i)+gamma)/((rA-rtol)*sin(gamma))));
        set(cone2,'Ydata',(1+sin(phidata(i)-gamma).*xrange./((rA-rtol)*sin(gamma)))./(cos(phidata(i)-gamma)/((rA-rtol)*sin(gamma))));
    end
    set(port,'Position',[rA*cos(phidata(i))-rA/4, rA*sin(phidata(i))-rA/4,rA/2,rA/2])
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