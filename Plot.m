%% simulate
%Plot Path

x = linspace(-3,10);
phi = params.phi; gamma = params.gamma; rp = params.rp; rtol = params.rtol;
y1 = (sin(phi+gamma).*x./((rp-rtol)*sin(gamma))-1)./(cos(phi+gamma)/((rp-rtol)*sin(gamma)));
y2 = (1+ sin(phi-gamma).*x./((rp-rtol)*sin(gamma)))./(cos(phi-gamma)/((rp-rtol)*sin(gamma)));
y1end = (sin(xtot(end,end)+gamma).*x./((rp-rtol)*sin(gamma))-1)./(cos(xtot(end,end)+gamma)/((rp-rtol)*sin(gamma)));
y2end = (1+ sin(xtot(end,end)-gamma).*x./((rp-rtol)*sin(gamma)))./(cos(xtot(end,end)-gamma)/((rp-rtol)*sin(gamma)));
figure(1)
plot(x,y1,'k','linewidth',2)
hold all
plot(x,y2,'k','linewidth',2)
plot(x,y1end,'k--','linewidth',2)
plot(x,y2end,'k--','linewidth',2)
rectangle('Position',[-rp,-rp,2*rp,2*rp],'Curvature',[1,1])
xlim([-1,6])
ylim([-1,6])
axis('square')
rectangle('Position',[rp*cos(phi)-rp/20, rp*sin(phi)-rp/20,rp/10,rp/10],'Curvature',[1,1])

% eta = params.eta; etaR = params.etaR; beta1 = params.beta1; beta2 = params.beta2; 
% beta3 = params.beta3; beta4 = params.beta4; beta5 = params.beta5;
N = params.N; Nc = params.Nc;
figure(1)
viscircles([xtot(1,end) xtot(2,end)],params.rs,'linewidth',1,'color','g');
viscircles([0 0],rp,'linewidth',1,'color','k');
hold all
viscircles([rp*cos(phi) rp*sin(phi)],rp/20,'color','k');
hold all
plot(xtot(1,:),xtot(2,:),'k')
plot((rp+rs)*cos(phi),(rp+rs)*sin(phi),'ro')
xlabel('X'); ylabel('Y');
%%
%Plot velocities and bounds
maxterm = abs(xtot(1,:)-(rp+rs)*cos(phi))+abs(xtot(2,:)-(rp+rs)*sin(phi));
%maxterm = abs(xtot(1,:)-rp*cos(phi))+abs(xtot(2,:)-rp*sin(phi))-(rp+rs)*(cos(phi)+sin(phi));

figure(9)
hold all
subplot(2,1,1)
hold all
plot(cos(xtot(3,:)+phi).*xtot(4,:)+sin(xtot(3,:)+phi).*xtot(5,:),'k')
plot(maxterm+beta1,'m--','linewidth',2)
legend('Normal Velocity','Normal Vel Bounds','location','southeast')
title(['Final Normal Velocity is ', num2str(cos(xtot(3,end)+phi).*xtot(4,end)+sin(xtot(3,end)+phi).*xtot(5,end)),' m/s'])
ylabel('speed [m/s]')
plot(-(maxterm+beta2),'m--','linewidth',2)
subplot(2,1,2)
hold all
plot(-sin(xtot(3,:)+phi).*xtot(4,:)+cos(xtot(3,:)+phi).*xtot(5,:),'k')
%plot(-sin(phi).*xtot(4,:)+cos(phi).*xtot(5,:),'k')
%plot(-sin(xtot(3,:)).*xtot(4,:)+cos(xtot(3,:)).*xtot(5,:),'k')
plot(maxterm+beta3,'g--','linewidth',2)
legend('Tangential Velocity','Tangential Vel Bounds','location','southeast')
title(['Final Tangential Velocity is ', num2str(-sin(xtot(3,end)+phi).*xtot(4,end)+cos(xtot(3,end)+phi).*xtot(5,end)),' m/s'])
plot(-(maxterm+beta3),'g--','linewidth',2)
xlabel('time')
ylabel('speed [m/s]')

%%
figure(3)
%plot(utot(4,:),'r')
hold all
%plot(utot(5,:),'b')
%plot(utot(6,:),'k')
plot(utot(10,:),'m')
plot(utot(11,:),'g')
legend('S_{Vt}', 'S_\Omega','location','best')
xlabel('time')
ylabel('Slack Variables')
title('Slack Variable Values')

%%
figure(4)
plot(xtot(6,:),'k')
hold all
plot(maxterm+beta4,'b--','linewidth',2)
hold all
legend('Rotational Speed','Rotational Speed Bounds')
plot(-(maxterm+beta5),'b--','linewidth',2)

xl = xlabel('time');
yl = ylabel(['$$\bf{\dot{\theta}}$$', '[rad/s]']);
set(yl, 'Interpreter', 'latex');
tit = title(['Final Rotational Speed is ',num2str(xtot(6,end)),' rad/s']);
% xlim([62 67])
% ylim([-1 1])

%%
% figure(6)
% plot(abs(xtot(1,:)-rp*cos(phi))+abs(xtot(2,:)-rp*sin(phi))-(rp+rs)*(cos(phi)+sin(phi))+phi*1.8,'g','linewidth',2)
% hold all
% plot(-(abs(xtot(1,:)-rp*cos(phi))+abs(xtot(2,:)-rp*sin(phi))-(rp+rs)*(cos(phi)+sin(phi))-phi/10),'g--','linewidth',2)
% plot(xtot(3,:))
figure(5)
plot(xtot(3,:),'k')
xlabel('time')
ylabel('\theta-\phi [rad]')
title(['Final Angle Error is ',num2str(xtot(3,end)),' rad'])
%%
figure(6)
plot(utot(1,:),'linewidth',2)
hold all
plot(utot(2,:),'linewidth',2)
plot(utot(3,:),'linewidth',2)
legend('u_x','u_y','u_{zz}','location','best')
xlabel('time')
ylabel('Inputs')
title('Thrust Inputs')


