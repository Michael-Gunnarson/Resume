clear
clc
format long

% EAE 129 Winter 2021 final project Michael Gunnarson


% Learjet C-21
Vinf = 170; %ft/s
W = 13000; %lb 
b = 34.0; %ft, span
cgbar = 0.32; %no units, cneter of gravity
Ixx = 2.79E4;
Izz = 4.11E4;  %^ Slug/ft^2
Xu = -0.0589; %1/s
Zu = -.3816; %1/s
ZdeltaE = -7.8612; %ft/s^2
Madot = -.3024; %1/s
MdeltaE = -2.8786; %1/s^2
q = 34.3; %lb/ft^2 dynamic pressure
S = 230; %ft^2
cbar = 7.0; %ft, wing chord
atrim = 5.0; %deg
Iyy = 1.88E4; %slug/ft^2
Ixz = -3.60E2; %ft/s^2
Xa = 11.3335; %ft/s^2
Za = -103.4862; %ft/s^2
Ma = -1.9387; %1/s^2
Mq = -0.8164; %1/s
Mu = -0.0002; %1/(ft*s)


U1 = Vinf;
% if it does not appear in chart assume zero
g = 32.174; % https://en.wikipedia.org/wiki/Standard_gravity


% Transfer Functions
s = tf('s');
% not given so assumed zero:
Xtu = 0;
BigTheta1 = 0;
Zadot = 0;
Zq = 0;
Mtu = 0;
Mta = 0;
% tf
A = [s-Xu-Xtu, -Xa, g*cos(BigTheta1); -Zu, s*(U1-Zadot)-Za, -(Zq+U1)*s+g*sin(BigTheta1); -(Mu+Mtu), -(Madot*s + Ma + Mta), s^2-Mq*s];
XdeltaE = 0;
B = [XdeltaE; ZdeltaE; MdeltaE];

TF = A\B;
TF = minreal(TF);
TFu = TF(1);
TFa = TF(2);
TFtheta = TF(3);

TFu = zpk(TFu)             %print me! 3.1 + char eqn 3.2
TFa = zpk(TFa)
TFtheta = zpk(TFtheta)

%datau = 
damp(TFu)           %print me! part of 3.2
%dataa = 
damp(TFa)
%datatheta = 
damp(TFtheta)

%[freq,damp,poleLocation] = damp(TFa)
% commented out because matlab doesn't like me assigning damp(TF-) twice

figure
impulse(TFu)
title('u(t) from Elevator Impulse Input')
ylabel('Perturbation velocity (ft/s)')

figure
step(TFu)
title('u(t) from Elevator Step Input')
ylabel('Perturbation velocity (ft/s)')

figure
impulse(TFa)
title('α(t) from Elevator Impulse Input')
ylabel('Perturbed Angle of Attack (deg)')

figure
step(TFa)
title('α(t) from Elevator Step Input')
ylabel('Perturbed Angle of Attack (deg)')

figure
impulse(TFtheta)
title('θ(t) from Elevator Impulse Input')
ylabel('Perturbed Angle W.R.T horizontal (deg)')

figure
step(TFtheta)
title('θ(t) from Elevator Step Input')
ylabel('Perturbed Angle W.R.T horizontal (deg)')


% Reduced order analysis

% s^2 - (Mq + Za/U1 + Madot)*s + (Za*Mq/U1-Ma) = 0

ReducedOrderChar = s^2 - (Mq + Za/U1 + Madot)*s + (Za*Mq/U1-Ma)

[omeg, zeta, pole] = damp(TFu)
shortPeriodNatFreq = omeg(3)
shortPeriodDampingRatio = zeta(3)
shortPeriodPole = pole(3)
PhugoidNatFreq = omeg(1)
PhugoidDampingRatio = zeta(1)
PhugoidPole = pole(1)
%short period wn and zeta > phugoid, index 3 or 4 from running code once
%and checking

% ShortPeriodChar = (s-pole(3))^2
% % didn't work
% PhugoidChar = (s-pole(1))^2

ShortPeriodChar = (s^2+2*zeta(3)*omeg(3)*s+omeg(3)^2)

PhugoidChar = (s^2+2*zeta(1)*omeg(1)*s+omeg(1)^2)

p = [1, -(Mq + Za/U1 + Madot), (Za*Mq/U1-Ma)];
reducedRoots = roots(p);
wn = (Za*Mq/U1-Ma)^0.5;
squigg = -(Mq + Za/U1 + Madot)/(2*wn);
% compare me to model char!!!



TF32 = (omeg(3)^2)/(s^2+2*zeta(3)*omeg(3)*s+omeg(3)^2)
TF342 = (wn^2)/(s^2+2*squigg*wn*s+wn^2)

figure
impulse(TF32)
title('Short Period Model Impulse Response')

figure
impulse(TF342)
title('Short Period Reduced Order Analysis Impulse Response')
