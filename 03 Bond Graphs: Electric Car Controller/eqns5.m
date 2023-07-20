function[ds, ext] = eqns5(t,s)
% CONSTANTS
global Rw Lw Tm M bt R Gr Cr g Cd rho Af T1 Kp Ki swap dt
% STATE VARIABLES
p3 = s(1);
p9 = s(2);
d = s(3);
dref = s(4);
if swap == "Part 1"
   d_dot = 0;
   dref_dot = 0;
   vref = 0;
% u in step of 100 volts for test:
if t >= 0 && t < T1
   uin = 0;
elseif t >= T1
   uin = 100; % volts, test case
end
elseif swap == "Part 2"
if t >= 0 && t < T1
   vref = 0;
elseif t >= T1
   vref = 1;
end
uin = Kp*(vref-p9/M) + Ki*(dref - d);
elseif swap == "Part 3"
vref = LA92Oracle(t);
uin = Kp*(vref-p9/M) + Ki*(dref - d);
end
%equations
f9 = p9/M;
% f9 = f10 = f11
d_dot = f9;
dref_dot = vref;
n = 0.0001;
sgn = @(v) v/(abs(v)+n);
e10 = 0.5*rho*Af*Cd*f9*abs(f9);
e11 = M*g*Cr*sgn(f9);
p3_dot = uin - (Rw/Lw)*p3 - (Tm*Gr/(R*M))*p9;
p9_dot = (Gr/R)*(Tm/Lw*p3 - bt/M*Gr/R*p9) - e10 - e11;
ds = [p3_dot;p9_dot;d_dot;dref_dot;uin];


%Pin = uin*p3/Lw;
%Pout = (p9/M)*p9/dt;
iin = p3/Lw;
vout = p9/M;

ext = [uin,p9_dot,iin,vout];
