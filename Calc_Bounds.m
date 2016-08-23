% Plot of max vtangent and min vnormal expected for successful grasping
collisionMomentum2D
%%
vt = mag_Tz.*sind(theta);
vn = mag_Jx.*cosd(theta);
figure
y = linspace(0,2);
x = linspace(0,max(vt)+max(vt)/4);
plot(x,vn)
hold all
plot(vt,y)
plot(
xlim([0,max(vt)+max(vt)/4])
ylim([0,2])
xlabel('Tangential velocity')
ylabel('Normal velocity')
legend('Minimium normal velocity','Maximum tangential velocity')
