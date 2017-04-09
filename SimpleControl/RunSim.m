ParamsS.rs = .15;
ParamsS.umax = .4/18.08;
ParamsS.tmas = 1;
ParamsS.xinit = [1 1 .05 -.1 pi/2 0]';
ParamsS.vf = .02;
ParamsT.rt = .15;
ParamsT.omega = -200*pi/180;
ParamsT.nu0 = -pi/2;
Ts1MPC = .5;
Ts1Simple = .01;
Ts2 = .01;
N = 20;
tolMPC = .001;
tolSimple = .05;

animate = 1;
%%
Tset = CalcTSet(ParamsS,Ts1MPC);

%%
[times, states, inputs] = SimMPC(ParamsS,ParamsT,Ts1MPC,Ts2,N,tolMPC,Tset,1);
%times(end)
%%
%AnimatePath;
[times,states,inputs,xdesr] = SimSimpleThrust(ParamsS,ParamsT,Ts1Simple,Ts2,tolSimple,1);
% figure
% idx = find(diff(times)>Ts1Simple+.01,1);
% times2 = 1:length(xdesr);
% times2 = times2.*Ts2;
% plot(times2,xdesr)
% hold all
% plot(t,xs)
figure(1)
hold all
plot(times,states(:,4))
%%
figure
plot(times,states(:,3))
%%
figure
plot(times,states(:,1))
hold all
plot(times,states(:,2))
plot(times,states(:,3))
