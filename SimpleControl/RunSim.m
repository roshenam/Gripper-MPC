ParamsS.rs = .15;
ParamsS.umax = .4/18;
ParamsS.tmas = 1;
ParamsS.xinit = [1 1 .1 .2 0 0]';
ParamsS.vf = .2;
ParamsT.rt = .15;
ParamsT.omega = -200*pi/180;
ParamsT.nu0 = -pi/2;
Ts1 = .5;
Ts2 = .01;
N = 20;
tol = .001;
animate = 1;
%%
Tset = CalcTSet(ParamsS,Ts1);

%%
[times, states] = SymMPC(ParamsS,ParamsT,Ts1,Ts2,N,tol,Tset,2);


%%
figure
plot(times,states(:,3))