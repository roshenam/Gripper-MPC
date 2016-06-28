N = [5 10 20 30];
Nc = [1 4 7 10];
Q = [10 10^2 10^3 10^4];
R = Q;

x0 =  5; y0 = 5; vx0 = 0; vy0 = 2;
init = [x0 y0 vx0 vy0];
phi = pi/4; omega = 5*pi/180;
params.phi = phi;
params.rp = .2; params.rtol = .2+.01;  params.rs = .2; params.gamma = pi/20;
params.Umax = .2;
params.Ts = .2;
params.Nc = 5;
results = cell(4,5);
for i=1:4
    params.N = N(i);
    [xtot, utot, cost, time] = MPC_Rotate_slack(init,params,phi,omega);
    results{i,1} = params;
    results{i,2} = xtot;
    results{i,3} = utot;
    results{i,4} = cost;
    results{i,5} = time;
end
%%
vals = zeros(4,7);
for i=1:4
    vals(i,1) = results{i,1}.N;
    vals(i,2) = sum(results{i,3}(3,:))+sum(results{i,3}(4,:));
    vals(i,3) = mean(results{i,5});
    vals(i,4) = std(results{i,5});
    vals(i,5) = sum(sqrt(results{i,3}(1,:).^2 + results{i,3}(2,:).^2));
    vals(i,6) = length(results{i,2});
    vals(i,7) = (results{i,2}(end,end) - atan2(results{i,2}(2,end),results{i,2}(1,end)))*180/pi;
end
outputs = struct('N',vals(:,1),'Sum_slack',vals(:,2),'Avg_comp_time',vals(:,3),'Std_Dev_comp_time',vals(:,4),'Total_Fuel',vals(:,5),'Total_Steps',vals(:,6),'Final_angle_error',vals(:,7));

