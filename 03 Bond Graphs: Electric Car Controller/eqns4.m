function[ds, ext] = eqns5(t,s)

% CONSTANTS

global Rw Rw Tm M bt R Gr Cr g Cd rho Af 

% STATE VARIABLES

p3 = s(1);
p9 = s(2);

%TIME STAMPS FOR u(t) AKA uin




%EQUATIONS

n = 0.001

p3_dot = uin - (Rw/Lw)*p3 - (Tm*Gr/(R*M))*p9;
p9_dot = (Gr/R)*((Tm/Lw)*p3 - (Gr*bt/(R*M))*p9) - (1/2)*rho*Af*Cd*(p9/M)*abs(p9/M) - Mg*Cr*(p9/M)/(p9/M + n)
