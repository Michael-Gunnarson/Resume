% EAE 126
% Spring 2022
% Project 2

% code written by:
% Michael Gunnarson

format shortG
clear
clc

% note: where a is considered the cylinder radius in calculations, R is
% the notation used in functions

%% Problem 1
% plot geometry of airfoils with various input parameters
% Joukowsky transformation

a = 1; % unit circle with radius of 1

Problem1(a)

%% Problem 2
% Find X component of surface velocity of airfoil
% plot surface pressure distrubution numerically (analytically is ec)

a = 1;
Uinf = 1;
tau = 0.1;

problem2(a,Uinf,tau)


%% Functions

function Problem1(a)
% plot 5 different joukowsky transformations of various circles

% plot first four graphs
e =  [0,0,-0.1,-0.1];
mu = [0,0.1,0,0.1];
title = ["Flat Plate","Circular Arc","Symmetrical Airfoil", "Cambered Airfoil"];

for i = 1:4
    plotAirfoil(e(i),mu(i),a,title(i))
end

% plot 5
e = 0;
mu = 0;
tau = 0.1;
b = b_solv2(a,tau);
title = "Ellipse";

plotAirfoil(e,mu,a,title,b)
end

function problem2(a,Uinf,tau)
% plotting U and Cp for each given case in problem 2

alpha = [0,5,5,5,0,0,5]*pi/180; % radians
k = ['n','n','y','y','n','n','y'];
e = [0,0,0,0,0,-0.1,-0.1];
mu = [0,0,0,0,0.1,0,0.1];

for i = 1:7
    if i >=1 && i <= 3
        b = b_solv2(a,tau);
        AirfoilSurfacePressure(Uinf,e(i),mu(i),a,alpha(i),k(i),b);
    else
    AirfoilSurfacePressure(Uinf,e(i),mu(i),a,alpha(i),k(i));
    end
end
end

function AirfoilSurfacePressure(Uinf,e,mu,a,alpha,k,b)
% find surface pressure and velocity of airfoil, plot

if ~exist('b','var')% check if b exists
    b = b_solv(e,a,mu); % solve b
end


% resolution
%r = 30;
theta = 100;
r = a;
R = a;

% check kutta condition
if k == 'n'
    gamma = gammaFunc(Uinf,R);
    t = "";
    len = 2*pi;
elseif k == 'y'
    gamma = gammaKutta(R,Uinf,alpha,mu);
    t = " with Kutta Condition";
    len = 2*pi-alpha; 
end

%r = linspace(R,rmax,r);
theta = linspace(0,len,theta);

% [r, theta] = meshgrid(r_data,theta_in);

% translate from polar to cartesian for plots (with offset)
x = r.*cos(theta) + e;
y = r.*sin(theta) + mu;


Vtheta = tangentialVelocity(Uinf,R,r,theta,gamma,alpha);

X = JX(x,y,b); % transforms x and y data
Y = JY(x,y,b);


U = surfaceVelocity(R,Vtheta,theta,X);

Cp = JoukowskyCp(U,Uinf);

% plot
figure

subplot(2,2,3)
plot(X(2:end),U)
title(sprintf("Surface Velocity at α = %g deg" + t,alpha*180/pi))
xlabel('X')
ylabel('U')

subplot(2,2,4)
plot(X(2:end),Cp)
title(sprintf("Surface Pressure Coefficient at α = %g deg" + t,alpha*180/pi))
xlabel('X')
ylabel('Cp')


subplot(2,2,1)
hold on
plot(x,y,'b')
%plot(x,y2,'b')
hold off
axis equal
title('Circular Cylinder')



subplot(2,2,2)
hold on
plot(X,Y,'b')
%plot(X,Y2,'b')
hold off
axis equal
title('Joukowsky Transformation')

end

function plotAirfoil(e,mu,a,t,b)
% plot resulting airfoil from data given in problem

if ~exist('b','var')% check if b exists
    b = b_solv(e,a,mu); % solve b
end

x = linspace(-a,a,500); % generate x data
[y,y2] = circXInput(x,a);

% offset:
x = x + e;
y = y + mu;
y2 = y2 + mu;

% plot circ

figure
sgtitle(sprintf(t+" ε = %g, µ = %g, b = %g",e,mu,b))

subplot(1,2,1)
hold on
plot(x,y,'b')
plot(x,y2,'b')
hold off
axis equal



X = JX(x,y,b); % transforms x and y data
Y = JY(x,y,b);
Y2 = JY(x,y2,b);

subplot(1,2,2)
hold on
plot(X,Y,'b')
plot(X,Y2,'b')
hold off
axis equal

end

function X = JX(x,y,b)
% Joukowsky transf. x
X = x.*(1+b.^2./(x.^2+y.^2));
end

function Y = JY(x,y,b)
% Joukowsky transf. y
Y = y.*(1-b.^2./(x.^2+y.^2));
end

function b = b_solv(e,a,mu)
% (technically e + or -)
b = e + (a^2-mu^2)^0.5;
end

function b = b_solv2(a,tau)
% ellipse with specified thickness
b = ((tau+1)/(1-tau))^0.5*a;
end

function [y,y2] = circXInput(x,a)
% using x as input, calc y data on circle
y = (a^2-x.^2).^0.5;
y2 = -(a^2-x.^2).^0.5;
end 

function vr = radialVelocity(Uinf,R,r,theta,alpha)
if ~exist('alpha','var')% check if variable exists
    alpha = 0;
end
vr = Uinf.*cos(theta-alpha).*(1-(R./r).^2);
end

function vtheta = tangentialVelocity(Uinf,R,r,theta,gamma,alpha)
if ~exist('alpha','var')% check if variable exists
    alpha = 0;
end
vtheta = -Uinf.*sin(theta-alpha).*(1+(R./r).^2) + gamma./(2*pi*r);
end

function gamma = gammaFunc(Uinf,R)
gamma = 4*pi*Uinf*R;
end

function gamma = gammaKutta(R,Uinf,alpha,mu)
gamma = -4*pi*R*Uinf*sin(alpha+asin(mu/R));
end

function U = surfaceVelocity(R,Vthet,thet,X)
% calculate the surface velocity at every point along surface

% initialize:
U = zeros(length(X)-2,1);

for i = 2:length(X)-1
    U(i) = R*Vthet(i)*(thet(i+1)-thet(i-1))/(X(i+1)-X(i-1));
end
end

function Cp = JoukowskyCp(U,Uinf)
Cp = -2*(U-Uinf)/Uinf;
end