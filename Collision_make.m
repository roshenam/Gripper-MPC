function [xi, yi] = Collision_make(params, c, nu, draw)

alphas = -2:1:6;
rs = params.rs; w = params.w;
rp = params.rp; 
thetas = alphas.*pi/4;
nus1 = thetas(1:5) + nu; 
x1s = c(1) + rs.*cos(nus1);
y1s = c(2) + rs.*sin(nus1);

beta1 = atan(rs/(rs+w));
rbig = sqrt(rs^2+(rs+w)^2);
x2 = c(1)+rbig*cos(pi+nu-beta1);
y2 = c(2)+rbig*sin(pi+nu-beta1);
x3 = c(1)+rbig*cos(pi+nu+beta1);
y3 = c(2)+rbig*sin(pi+nu+beta1);

x = [x1s x2 x3]';
y = [y1s y2 y3]';

thetas2 = 0:pi/8:2*pi;
xcirc = rp.*cos(thetas2);
ycirc = rp.*sin(thetas2);

if draw
figure(1)
fill(x,y,'b')
hold all
fill(xcirc,ycirc,'r')
xl = xlim;
yl = ylim;
l = min(xl(1),yl(1));
u = max(xl(2),yl(2));
xlim([l,u])
ylim([l,u])
axis('square')
end
[xi,yi] = polyxpoly(x,y,xcirc,ycirc);

