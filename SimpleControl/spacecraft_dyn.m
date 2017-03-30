function dydt = spacecraft_dyn(t,y,P)

idx = find(P.t>=t,1);
if isempty(idx)
    idx = length(P.t);
end
u = P.u(idx,:);
u=u';
%P.i = P.i+1;
dydt = zeros(4,1);
xdot = P.A*y+P.B*u;
dydt(1) = xdot(1);
dydt(2) = xdot(2);
dydt(3) = xdot(3);
dydt(4) = xdot(4);


