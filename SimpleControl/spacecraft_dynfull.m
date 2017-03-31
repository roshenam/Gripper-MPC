function dydt = spacecraft_dynfull(t,y,P)

y_des = [P.vxd*t; P.vyd*t; P.vxd; P.vyd; 0; 0];
if (P.tdis0<t && t<P.tdisf)
    u_dis = P.udis;
else
    u_dis = 0;
end
u = -P.K*(y(1:6)-y_des)-P.Ki*y(7:8);
if u(1) > P.umax
    u(1) = P.umax;
elseif u(1) < -P.umax
    u(1) = -P.umax;
end
if u(2) > P.umax
    u(2) = P.umax;
elseif u(2) < -P.umax
    u(2) = -P.umax;
end
% figure(3)
% hold all
% plot(t,u(1),'ro')
% hold all
% plot(t,u(2),'go')
u = u+[u_dis; 0; 0];
dydt = zeros(8,1);
xdot = P.A*y+P.B*u+P.R*y_des;
dydt(1) = xdot(1);
dydt(2) = xdot(2);
dydt(3) = xdot(3);
dydt(4) = xdot(4);
dydt(5) = xdot(5);
dydt(6) = xdot(6);
dydt(7) = xdot(7);
dydt(8) = xdot(8);



