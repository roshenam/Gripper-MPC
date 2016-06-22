clear all
syms eta phi beta x y z X Y Z
R2 = [1 0 0; 0 cos(phi) sin(phi); 0 -sin(phi) cos(phi)];
R1 = [cos(eta) sin(eta) 0; -sin(eta) cos(eta) 0; 0 0 1];
R = R2*R1;
vIJK = R*[X;Y;Z];
values = (x.^2 + y.^2).^.5 - z.*tan(beta);
values = subs(values, x, vIJK(1));
values = subs(values, y, vIJK(2));
values = subs(values, z, vIJK(3));
values = simplify(values);
%%
x = linspace(-3,3);
y = linspace(-3,3);
z = linspace(-1,3);
beta = 10*pi/180;
eta = 0; phi = -pi/4;
[X,Y,Z] = meshgrid(x,y,z);
values = ((Z.*sin(phi) + Y.*cos(eta)*cos(phi) - X.*cos(phi)*sin(eta)).^2 +...
    (X.*cos(eta) + Y.*sin(eta)).^2).^(1/2) - tan(beta).*(Z.*cos(phi) -...
    Y.*cos(eta)*sin(phi) + X.*sin(eta)*sin(phi));


idx = find(values<=0);
plot3(X(idx),Y(idx),Z(idx),'ro','markersize',1)
xlabel('X')
ylabel('Y')
zlabel('Z')
grid on
lims = axis;
hold all
plot3(lims(1:2),[0 0],[0 0],'k') 
plot3([0,0],lims(3:4),[0,0],'k')
plot3([0,0],[0,0],lims(5:6),'k')

[xs, ys, zs] = sphere(50);
surf(.2.*xs, .2.*ys, .2.*zs)
plot3(0,3,2.5,'ko','markersize',4)
view([120 20])

