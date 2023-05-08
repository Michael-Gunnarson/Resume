% EAE 126, Spring 2022
% Project 8
% Code written by Michael Gunnarson

close all
%format shortG
clear
clc

% Nonlinear Supersonic and Hypersonic Flow

%% Problem 1
% wedge at alpha = 0.  Find shock angle and Cp vs deflection angle for Minf
% = 2 4 5 8 10 inf
% compare with newtonian results

gamma = 1.4;

Ma_vec = [2,4,6,8,10,inf];
theta = [0,5,10,15,20]*pi/180; % deg

figure
hold on

for i = 1:length(Ma_vec)
    Ma = Ma_vec(i);
beta = findBeta(gamma,Ma,theta);
cp = findCp(beta,gamma,Ma);


plot(theta,beta)
end

ylabel('beta')
xlabel('theta (rads)')
title(sprintf('beta vs theta'))% for Ma = %g',Ma))
legend('Ma = 2','Ma = 4','Ma = 6','Ma = 8','Ma = 10','Ma = inf')
hold off


figure
hold on

for i = 1:length(Ma_vec)
    Ma = Ma_vec(i);
beta = findBeta(gamma,Ma,theta);
cp = findCp(beta,gamma,Ma);


plot(theta,cp)
end

ylabel('cp')
xlabel('theta (rads)')
title(sprintf('Cp vs theta'))% for Ma = %g',Ma))
legend('Ma = 2','Ma = 4','Ma = 6','Ma = 8','Ma = 10','Ma = inf')
hold off


%% Problem 2
% flat plate at alpha = 2deg.  Calculate lift and drag coefficients from
% shock expansion theory, Ma = 2 4 6 8 10 inf
% compare with Newtonian Theory

p1 = 100; % kPa
b = 1; % meter
c = b;
gamma = 1.4;
alpha = 2*pi/180; % rads

Ma_vec = [2,4,6,8,10,inf];
for i = 1:length(Ma_vec)
Ma = Ma_vec(i);

% fan
p0 = isentropicStagnation(p1,Ma);
w = omegaFunction(gamma,Ma);
dw = alpha;
w2 = w + dw;
Ma2 = omegaFunctionReverse(gamma,w2);
p3 = isentropicStagnation2(p0,Ma2);


% lower oblique shock
theta = alpha;
beta = findBeta(gamma,Ma,theta);
Man1 = sin(beta)*Ma;
Man2 = normalMach2(Man1,gamma);
p2_over_p1 = pressureRelation(Ma,beta,gamma);
Ma2_lower = Man2/sin(beta-theta);
% find p2 with ratio
p2 = p1*p2_over_p1;
% compute lift and drag coefficients

CL(i) = liftCoeff(alpha,p2,p3,b,c,gamma,p0,Ma);
CD(i) = dragCoeff(alpha,p2,p3,b,c,gamma,p0,Ma);

%disp(i)
end

figure
plot(Ma_vec,CL)
title('Lift Coefficient vs Mach Number')
xlabel('Mach Number')
ylabel('CL')

figure
plot(Ma_vec,CD)
title('Drag Coefficient vs Mach Number')
xlabel('Mach Number')
ylabel('CD')


%% Functions

function beta = findBeta(gamma,Ma,theta)
beta = theta.*((gamma+1)/4+(((gamma+1)/4)^2+1./(Ma^2*theta.^2)).^0.5);
end

function cp = findCp(beta,gamma,Ma)
cp = 4/(gamma+1)*((sin(beta)).^2-1./Ma.^2);
end

function p02_over_p01 = stagnationPressure(k,Ma,beta)
p02_over_p01 = (((k+1)*Ma^2*(sin(beta)).^2)/(2+(k-1)*Ma^2*(sin(beta)).^2)).^(k/(k-1))*((k+1)/(2*k*Ma^2*(sin(beta)).^2-(k-1))).^(1/(k-1));
end

function p02 = stag2(p01,k,Ma,beta)
p02 = p01*stagnationPressure(k,Ma,beta);
end

function p01 = stag1(p02,k,Ma,beta)
p01 = p02/stagnationPressure(k,Ma,beta);
end

function p0 = isentropicStagnation(p1,Ma)
p0 = p1*(1+0.2*(Ma).^2).^3.5;
end

function p1 = isentropicStagnation2(p0,Ma)
p1 = p0/(1+0.2*(Ma).^2).^3.5;
end

function w = omegaFunction(gamma,Ma)
K = (gamma+1)/(gamma-1);
w = K^0.5*atan((Ma^2-1)/K)^0.5-atan(Ma^2-1)^0.5;
end

function Ma2 = omegaFunctionReverse(gamma,w)

res = 0.5;

% initial 2 mach # guesses for Newton's method
Ma2 = 1.5;
Ma_new = 1.7;



g = @(Ma) omegaFunction(gamma,Ma)-w; % desired input of Ma results in g=0

while res > 1E-6

    % cycle
    Ma = Ma2;
    Ma2 = Ma_new;

g1 = g(Ma);
g2 = g(Ma2);

dg = (g2-g1)/(Ma2-Ma);
Ma_new = Ma - g1/dg;

res = abs(Ma_new-Ma2);
end
end

function L = lift(alpha,p2,p3,b,c)
L = cos(alpha)*force(p2,p3,b,c);
end

function D = drag(alpha,p2,p3,b,c)
D = sin(alpha)*force(p2,p3,b,c);
end

function F = force(p2,p3,b,c)
F = (p2-p3)*b*c;
end

function CL = liftCoeff(alpha,p2,p3,b,c,gamma,pinf,Ma)
CL = lift(alpha,p2,p3,b,c)/(0.5*gamma*pinf*Ma^2*b*c);
end

function CD = dragCoeff(alpha,p2,p3,b,c,gamma,pinf,Ma)
CD = drag(alpha,p2,p3,b,c)/(0.5*gamma*pinf*Ma^2*b*c);
end

function Man2 = normalMach2(Man1,k)
Man2 = (((k-1)*Man1^2+2)/(2*k*Man1^2-(k-1)))^0.5;
end

function p2_over_p1 = pressureRelation(Ma,beta,k)
p2_over_p1 = 1/(k+1)*(2*k*Ma^2*(sin(beta))^2-(k-1));
end