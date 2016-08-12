function [out] = Inv_Set_Sim(params,loop)
% Defining parameters
x0vec = params.x0vec;
omega = params.omega;
rp = params.rp; rs = params.rs;
Ts = params.Ts;

% Simulation
Nsim = params.Nsim;
data = loop.simulate(x0vec(1:6), Nsim);
[~,n] = size(data.X);

if n==1
    out.X = NaN; out.U = NaN;
    out.XnonI = NaN; out.UnonI = NaN;
    out.final = NaN;
    return
end

xtotnonI = data.X; % Save data in non-inertial frame
utotnonI = data.U; % Save thrust data in non-inertial frame
xtot = zeros(7,Nsim+1);
utot = zeros(3,Nsim);
mult = 0:Nsim;

phi = params.phi;
xtot(7,:) = phi+Ts*omega.*mult;
xtot(1,:) = xtotnonI(1,:).*cos(xtot(7,:))-xtotnonI(2,:).*sin(xtot(7,:));
xtot(2,:) = xtotnonI(1,:).*sin(xtot(7,:))+xtotnonI(2,:).*cos(xtot(7,:));
xtot(3,:) = xtotnonI(3,:) + xtot(7,:);
xtot(4,:) = xtotnonI(4,:).*cos(xtot(7,:))-xtotnonI(5,:).*sin(xtot(7,:));
xtot(5,:) = xtotnonI(4,:).*sin(xtot(7,:))+xtotnonI(5,:).*cos(xtot(7,:));
xtot(6,:) = xtotnonI(6,:) + omega;

utot(1,:) = utotnonI(1,:).*cos(xtot(7,1:end-1))-utotnonI(2,:).*sin(xtot(7,1:end-1));
utot(2,:) = utotnonI(1,:).*sin(xtot(7,1:end-1))+utotnonI(2,:).*cos(xtot(7,1:end-1));
utot(3,:) = utotnonI(3,:);

% Determining incoming conditions upon impact
dist = sqrt(data.X(1,:).^2+data.X(2,:).^2);
idx = find(dist<=(rp+rs),1);
if isempty(idx)
    final = NaN;
else
    final.idx = idx;
    final.vn = data.X(4,idx-1);
    final.vt = data.X(5,idx-1);
    final.vmag = norm([data.X(4,idx-1) data.X(5,idx-1)]);
    final.nu = data.X(3,idx-1)*180/pi;
    final.angattack = atan2(data.X(5,idx-1),data.X(4,idx-1));
    final.nudot = data.X(6,idx-1)*180/pi;
    final.thrust = sum(sqrt(utot(1,:).^2+utot(2,:).^2));
end
out.X = xtot; out.U = utot; out.XnonI = xtotnonI; out.UnonI = utotnonI;
out.final = final;

