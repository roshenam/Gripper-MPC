collisionMomentum2D

%%
vt = mag_Tz.*sind(theta);
vn = mag_Jx.*cosd(theta);
%%
x = linspace(0,max(vt)+max(vt)/3);
y = linspace(0,2.2);
plot(x,vn)
hold all
plot(vt,y)
plot(x,2.*ones(size(x)))
ylim([0,2.2])
