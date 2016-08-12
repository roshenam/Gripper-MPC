% Run simulations across variations of Q and R and N with varying x, y, xdot,
% ydot to check region of attraction, cost variations, and final conditions
params.phi = pi/4;
params.rp = .2; params.rtol = .2+.01;  params.rs = .2; params.gamma = 20*pi/180;
params.Umax = .2;
params.Tmax = 1;
params.Ts = .2;
params.betaHIGH = -1.5;
params.eta = 1;
params.tantol = .4;
params.Nsim = 70;

weights = [10^2 10^4; 10^3 10^3; 5*10^3 5*10^2; 10^4 10^2];

counter = 0;
for omega = [2*pi/180 5*pi/180]
    for i=1:4
        for n = [10 20 30]
            counter = counter + 1;
            params.N = n;
            params.omega = omega;
            params.Qval = weights(i,1);
            params.Rval = weights(i,2);
            loop = Inv_Set_Calc(params);
            tests(counter).loop = loop;
            tests(counter).params = params;
        end
    end
end

%%
xs = linspace(0,3,5); ys = [-.75, -.5, -.25, .25, .5, .75]; vxs = [-.9, -.5, -.2, 0, .2, .5, .9];
vys = vxs;
counter = 0;
failures = [];
for i = 1:24
    for x = xs
        for y = ys
            for vx = vxs
                for vy = vys
                    counter = counter + 1;
                    loop = tests(i).loop;
                    params = tests(i).params;
                    x0 = x; y0 = y*tan(params.gamma); theta0 = -pi/4;
                    vx0 = vx; vy0 = vy; thetadot0 = 5*pi/180;
                    Rmat = [cos(params.phi) -sin(params.phi); sin(params.phi) cos(params.phi)];
                    rI = Rmat*[x0;y0];
                    vI = Rmat*[vx0;vy0];
                    init = [rI(1) rI(2) theta0 vI(1) vI(2) thetadot0];
                    nu0 = theta0-params.phi; nu0dot = thetadot0-params.omega;
                    x0vec = [x0; y0; nu0; vx0; vy0; nu0dot];
                    params.x0vec = x0vec;
                    params.init = init;
                    DATA(counter).failure = 0;
                    try
                    [data] = Inv_Set_Sim(params,loop);
                    catch
                        DATA(counter).failure = 1;
                        DATA(counter).data = NaN;
                        DATA(counter).params = params;
                        DATA(counter).loop = loop
                        continue
                    end
                    DATA(counter).data = data;
                    DATA(counter).params = params;
                    DATA(counter).loop = loop;
                    disp(['Finished simulation ',num2str(counter)])
                end
            end
        end
    end
end


% For each test, save params structure and data structure
for i=1:4
    x0 =  states(i,1); y0 = states(i,2); theta0 = -pi/4; vx0 = states(i,3); vy0 = states(i,4); thetadot0 = 5*pi/180;
    Rmat = [cos(phi) -sin(phi); sin(phi) cos(phi)];
    rI = Rmat*[x0;y0];
    vI = Rmat*[vx0;vy0];
    init = [rI(1) rI(2) theta0 vI(1) vI(2) thetadot0];
    x0vec = [x0; y0; nu0; vx0; vy0; nu0dot];
    params.x0vec = x0vec;
    [data] = Inv_Set_Sim(params);
end
