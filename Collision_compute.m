function [vnorm, angattack, offset, collend] = Collision_compute(params,data)
for i=1:length(data.X)
    [xi, yi] = Collision_make(params, data.X(1:2,i), data.X(3,i), 0);
    if isempty(xi)
    else
        collend = i;
        collpre = i-1;
        break
    end
end
if isempty(xi)
    disp('Collision never occured')
end

pospre = data.X(1:2,collpre);
vend = data.X(4:5,collend);
vpre = data.X(4:5,collpre);
nuend = data.X(3,collend);
nupre = data.X(3,collpre);
vavg = [mean([vend(1);vpre(1)]); mean([vend(2);vpre(2)])];
nuavg = mean([nuend(1);nupre(1)]);
R = [cos(nuavg) sin(nuavg); -sin(nuavg) cos(nuavg)];
vrot = -R*vavg;
vnorm = norm(vrot);
angattack = atan2(vrot(2),vrot(1))*180/pi;

time = 0:params.Ts/10:params.Ts;
xinterp = pospre(1)+time.*vavg(1);
yinterp = pospre(2)+time.*vavg(2);
nuinterp = nupre + time.*mean([data.X(6,collpre),data.X(6,collend)]);

for i=1:length(time)
    [xi, yi] = Collision_make(params, [xinterp(i); yinterp(i)], nuinterp(i), 0);
    if isempty(xi)
    else
        break
    end
end
% 
% Collision_make(params, [data.X(1,collpre); data.X(2,collpre)],data.X(3,collpre), 1);
% title('pre collision')
 Collision_make(params, [xinterp(i); yinterp(i)],nuinterp(i), 1);
 title('fine collision')
% Collision_make(params, [data.X(1,collend); data.X(2,collend)],data.X(3,collend), 1);
% title('should be same as first')


collIfine = i;
Rcoll = [cos(nuinterp(collIfine)) sin(nuinterp(collIfine)); -sin(nuinterp(collIfine)) cos(nuinterp(collIfine))];
collpos = Rcoll*[mean(xi); mean(yi)];
centerlocx = xinterp(collIfine)-(params.rs+params.w)*cos(nuinterp(collIfine));
centerlocy = yinterp(collIfine)-(params.rs+params.w)*sin(nuinterp(collIfine));
centerpos = Rcoll*[centerlocx; centerlocy];
offset = abs(collpos(2)-centerpos(2));

