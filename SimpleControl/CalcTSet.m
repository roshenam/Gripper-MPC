function Tset = CalcTSet(ParamsS,Ts1)

A = zeros(4,4);
A(1,3) = 1; A(2,4) = 1;
B = zeros(4,2);
B(3,1) = 1; B(4,2) = 1;
C = zeros(2,4);
C(1,1) = 1; C(2,2) = 1;
D = zeros(2,2);
sysC = ss(A,B,C,D);
%Ts = .5;
sysD = c2d(sysC,Ts1);
Ad = sysD.A; Bd = sysD.B;
system = LTISystem('A',Ad,'B',Bd,'Ts',Ts1);
%umax = .4/18;
umax = ParamsS.umax;
Q = eye(4); Q(3,3) = 100; Q(4,4) = 100;
R = 10^2.*eye(2);
system.x.min = [-2; -2; -1; -1];
system.x.max = [2; 2; 1; 1];
system.u.min = [-umax; -umax];
system.u.max = [umax; umax];
system.x.penalty = QuadFunction(Q);
system.u.penalty = QuadFunction(R);
Pn = system.LQRPenalty;
system.x.with('terminalPenalty');
system.x.terminalPenalty = Pn;
%Tset = system.LQRSet;
Tset = system.invariantSet();
