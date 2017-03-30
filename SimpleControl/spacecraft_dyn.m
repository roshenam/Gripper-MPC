function dydt = spacecraft_dyn(t,y,P)

idx = find(P.t>=t,1);
if isempty(idx)
    idx = length(P.t);
end
unom = P.u(idx,:);
unom=unom';
y_des = P.state_des(idx,:);
%y_noise = randn(4,1)./10;
u_fric = -.01.*[y(3); y(4)]./sqrt(y(3)^2+y(4)^2);
uLQR = -P.K*(y-y_des');
utot = unom+uLQR;
if utot(1) > P.umax
    utot(1) = P.umax;
elseif utot(1) < -P.umax
    utot(1) = -P.umax;
end

if utot(2) > P.umax
    utot(2) = P.umax;
elseif utot(2) < -P.umax
    utot(2) = -P.umax;
end

%figure(20)
%plot(idx,uLQR(1),'ro')
%hold all
%plot(idx,uLQR(2),'go')
%u = utot+u_fric;
%u = unom+u_fric;
%u = unom+u_fric;
%u = unom;
u = utot;
%P.i = P.i+1;
dydt = zeros(4,1);
%xdot = P.A*(y+y_noise)+P.B*u;
xdot = P.A*y+P.B*u;
dydt(1) = xdot(1);
dydt(2) = xdot(2);
dydt(3) = xdot(3);
dydt(4) = xdot(4);


