% code written by 
% Michael Gunnarson
% as a part of the UC Davis EAE 143 Psyche team Winter 2021:
% Jackie Arroyo Donjuan, Nickolas Loftus, Michael Gunnarson, Christian Lum.


% Orbital mechanics hand calculations use patched conic sections to acquire
% rough estimates of transfer from 16 Psyche to 216 Kleopatra.
% Assumes instantaneous impulse and parallel exit velocity with planet.
% Assumes that timing is chosen st all asteroids will be in ideal positions
% Uses sphere of influence assumption.  Assumes 2D plane


% This code is used to calculate how much delta v we would need to get from
% 16 Psyche to 216 Kleopatra, both of which are large metallic asteroids.
% More delta v would be required to perform data maneuvers around 216
% Kleopatra.

% updated December 2022

clear
format long
clc


%% Constants

% physical and orbital data from Francis J. Hale's Introduction to
% Spaceflight.  Psyche data taken from https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html#/?sstr=16
% 216 Kleopatra data taken from https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html#/?sstr=2000216
% 16 denotes 16 psyche, 216 denotes 216 kleopatra

% radius
rEarth = 6.378E6; % m
rMars = 71.40E6; % m
rSun = 696.0E6;
r16 = 226/2*1000; % diameter
r216 = 54.25*1000; % equivalent radius, m, https://www.sciencedirect.com/science/article/pii/S0019103510004355

% mu
muEarth = 3.986e14; % m^3/s^2
muMars = 4.297e13; % m^3/s^2
muSun = 1.327E20;
mu16 = 1.53*1000^3; % m^3/s^2 

% mass to mu
m216 = 4.64E18; % kg https://www.sciencedirect.com/science/article/pii/S0019103510004355
G = 6.674E-11; % m^3/kgs^2 from https://physics.nist.gov/cgi-bin/cuu/Value?bg
mu216 = G*m216;

% semimajor axis
aEarth = 1.000; % AU
aMars = 1.524; % AU
aBelt = 2.06; % fromhttps://en.wikipedia.org/wiki/Asteroid_belt
a16 = 2.924484673106921; % AU
a216 = 2.793277353255009; % AU

% eccentricity
eEarth = 0.0167;
eMars = 0.0934;
e16 = 0.133926624302819;
e216 = .2507297431057924;

% orbital velocity
velEarth = 29790; % m/s
velMars = 24140; % m/s

% Asteroid Orbits
ra16 = 3.316151033201463; % AU
rp16 = 2.532818313012378; % AU
ra216 = 3.493635066459865; % AU
rp216 = 2.092919640050153; % AU

%% Solar-Centric Transfer

% unit conversions
a16 = AUtoM(a16);
ra16 = AUtoM(ra16);
rp16 = AUtoM(rp16);

a216 = AUtoM(a216);
ra216 = AUtoM(ra216);
rp216 = AUtoM(rp216);


[dv2,dv3,tTsfr,va16,va216,vp16,vp216,vaTsfr] = transfer(a216,ra216,rp216,e216,a16,ra16,rp16,e16,rSun,muSun);


%% 16 Psyche Escape
vInf = dv2;
rOrbitD = 191E3; % m
% from Psyche Science Operations Concept: Maximize Reuse to Minimize Risk

[dv1,tEsc] = escape(mu16,muSun,rOrbitD,vInf,r16,a16);


%% 216 Kleopatra Capture

vInf = abs(vaTsfr-va216);
[dv4,tCap,rp,vHyperAtrp] = capture(mu216,muSun,a216,r216,vInf);


%% The Verdict

t1Hr = StoH(tEsc);
t2Hr = StoH(tTsfr);
t3Hr = StoH(tCap);

t1D = StoD(tEsc);
t2D = StoD(tTsfr);
t3D = StoD(tCap);

tTotalS = tEsc + tTsfr + tCap;
tTotalYr = StoYr(tTotalS);

deltaV = dv1+dv4; % dv2 and 3 are left out because there is no firing required in the hohmann transfer

                      


sprintf('The escape phase requires %Gm/s and takes %G hr',dv1,t1Hr)
sprintf('The Hohmann transfer phase requires %Gm/s and takes %G days',0,t2D)
sprintf('The capture phase requires %Gm/s and takes %G hr',dv4,t3Hr)
sprintf('The entire mission (escape to capture) requires %Gm/s and takes %G years',deltaV,tTotalYr)
%sprintf('More ΔV is required to utilize payload in a similar fashion to Psyche Mission, but these results are quite heartening')
sprintf('We find that the required ΔV for full capture is larger than the budget')
sprintf('Budget is 783.9 m/s')


%% Functions

% E = -mu^2/(2H^2) % for circles only (eccentricity = 0)
                   % with thanks to Marti Sarigul-Klijn for the formula

function Vr = Oberth(mu,r,Vinf)
Vr = ((escapeVelocity(mu,r))^2+Vinf^2)^0.5;
end

function Vesc = escapeVelocity(mu,r)
Vesc = (2*mu/r)^0.5;
end

function Vcirc = circularVelocity(mu,r)
Vcirc = (mu/r)^0.5;
end

function E = energy1(mu,a)
E = -mu/(2*a);
end

function E = energy2(V,mu,r)
E = V^2/2-mu/r;
end

function Einf = energyInfinity(Vinf)
% where Vinf is the velocity at r = infinity
Einf = Vinf^2/2;
end


function p = semilatus1(a,epsilon) 
% semilatus rectum, distance above foci 
p = a*(1-epsilon^2);
end

function p = semilatus2(H,mu) 
% semilatus rectum, distance above foci 
p = H^2/mu;
end

function e  = eccentricity1(ra,rp)
e = (ra-rp)/(ra+rp);
end

function e  = eccentricity2(E,H,mu)
e = (1+(2*E*H^2/mu^2))^0.5;
end

function H = angularMomentum4(va,ra)
% works at both apogee and perogee
H = va*ra;
end


function a = majorAxisHalf1(ra,rp) % where 2a is the major axis
a = (ra+rp)/2;
end

function a = majorAxisHalf2(p,e) % where 2a is the major axis
a = p/(1-e^2);
end

function a = majorAxisHalf3(mu,E) % where 2a is the major axis
a = -mu/(2*E);
end

function rp = radiusPeriapsis(a,e) % smallest radius in orbit
rp = a*(1-e);
end

function V = velocity(E,mu,r) % general radius
V = (2*(E+mu/r))^0.5;
end


function x = parax(rho,theta)
x = rho.*cos(theta);
end

function y = paray(rho,theta)
y = rho.*sin(theta);
end

function [x,y] = para(rho,theta)
x = rho.*cos(theta);
y = rho.*sin(theta);
end

function [x,y] = paraHyperbola(a,p)
x = linspace(a, 2*a,1000);
b = bGeo(a,p);
c = (a^2+b^2)^0.5;
y = b*(x.^2./a^2-1).^0.5;
x = x+c;
end


function m = AUtoM(AU)
% from https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html#/?sstr=16&view=OPDA
% convert astronomical units to meters
m = AU*149597870.70*1000;
end

function r = SOI(m2,m1,a)
% where m2 is the smaller mass, m1 is the larger mass, and a is the
% semi-major axis of m2 revolving m1
r = a*(m2/m1)^0.4;
end

function b = bGeo(a,p)
b = (abs(a)*p)^0.5;
end

function t = TOF_ellipse(a,mu)
t = 2*pi*(a^3/mu)^0.5;
end

function nu = trueAnomoly(e,p,r)
nu = acos(1/e*(p/r-1));
end

function tof = TOF_hyp(nu,mu,e,a)
% to find the hyperbolic time of flight from periapsis to point A
% thus nu should be the true anomoly at point A

coshf = (e+cos(nu))/(1+e*cos(nu));
F = acosh(coshf);
M = e*sinh(F)-F;
tof = M*((-a)^3/mu)^0.5;

end

function hr = StoH(s)
hr = s/(60*60);
end

function hr = StoD(s)
hr = s/(60*60*24);
end

function hr = StoYr(s)
hr = s/(60*60*24*365);
end

function v = derivedV(mu,a,r)
% derived formula for v(mu,a,r) in order to minimize truncation error
v = (mu*(2*a-r)/(a*r))^0.5;
end

function b = impactParameter(vInf,mu,r)
vEsc = escapeVelocity(mu,r);
b = r*(1+vEsc^2/vInf^2)^0.5;
end

function [dv1,tEsc] = escape(mu,muSun,rOrbitD,Vinf,r16,a16)
% planet centric calculations

Vcirc = circularVelocity(mu,rOrbitD); % velocity of spacecraft orbit D
Vob = Oberth(mu,rOrbitD,Vinf); % velocity needed to escape 16 psyche
dv1 = abs(Vob-Vcirc); 

% orbit characteristics

E = energyInfinity(Vinf);
H = angularMomentum4(Vob,rOrbitD);
e = eccentricity2(E,H,mu);
p = semilatus2(H,mu);
a = majorAxisHalf2(p,e);

% TOF
soi = SOI(mu,muSun,a16);
nu = trueAnomoly(e,p,soi);
tEsc = TOF_hyp(nu,mu,e,a);

% graph

theta = linspace(0,2*pi,500);
[dx,dy] = para(rOrbitD,theta); % orbit D x and y
psy_y = paray(r16,theta); % psyche y data
[ex,ey] = paraHyperbola(a,p); % escape x and y

figure
title(sprintf('16 Psyche Escape from Orbit D Vinf = %g m/s, μ = %g',Vinf, mu))
hold on
plot(dx,dy,'b','Linewidth',2) % orbit D
for i = 1:100:r16 % plot psyche as cyan
    psy_x = parax(i,theta);
    plot(psy_x,psy_y,'c')
end
plot(ex,ey, 'r','Linewidth',2) % escape
xlabel('m')
ylabel('m')
set(gca, 'color', [0,0,0])
axis equal
hold off


end



function [dv1,dv2,tTsfr,va16,va216,vp16,vp216,vaTsfr] = transfer(a216,ra216,rp216,e216,a16,ra16,rp16,e16,rSun,muSun)
% for inner or outer hohman transfer orbit, from ellipse to ellipse
% sun centric calculations

E16 = energy1(muSun,a16);
E216 = energy1(muSun,a216);

check = E16 > E216;

switch check
    case 1 % fire opposite direction of planet
        
        % outer asteroid
        va16 = derivedV(muSun,a16,ra16);
        vp16 = derivedV(muSun,a16,rp16);
        p16 = semilatus1(a16,e16);
        
        
        % inner asteroid
        va216 = derivedV(muSun,a216,ra216);
        vp216 = derivedV(muSun,a216,rp216);
        p216 = semilatus1(a216,e216);

        % Hohman transfer orbit
        rpTsfr = rp16; % by nature of inner orbit geometry firing at apogee
        raTsfr = ra216; % by nature of inner orbit geometry firing at apogee
        % for some reason this only works with overlapping ellipses.  This
        % is the case we find ourselves in, but do not use this code
        % blindly later.  It appears swapping raTsfr and rpTsfr does the
        % trick but since I can't tell when it will overlap, I can't tell
        % wwhen it will work and when it won't
        
        
        eTsfr  = eccentricity1(raTsfr,rpTsfr);
        aTsfr = majorAxisHalf1(raTsfr,rpTsfr);
        pTsfr = semilatus1(aTsfr,eTsfr); 
        vaTsfr = derivedV(muSun,aTsfr,raTsfr);
        vpTsfr = derivedV(muSun,aTsfr,rpTsfr);

        tTsfr = TOF_ellipse(aTsfr,muSun)/2;
        t16 = TOF_ellipse(a16,muSun)/2;
      
 
        dv1 = abs(vp16-vpTsfr); % dv leaving psyche
        dv2 = abs(vaTsfr-va216); % dv to enter kleopatra


        theta = linspace(0,2*pi,500);

        bTsfr = bGeo(aTsfr,pTsfr);
        cTsfr = (aTsfr^2-bTsfr^2)^0.5;

        xTsfr = (aTsfr)*cos(theta)-cTsfr; % since a > b, focus (sun) is offset on the x axis
        yTsfr = (bTsfr)*sin(theta);

        b16 = bGeo(a16,p16);
        c16 = (a16^2-b16^2)^0.5;

        x16 = (a16)*cos(theta)-c16; % since a > b, focus (sun) is offset on the x axis
        y16 = (b16)*sin(theta);
        
        b216 = bGeo(a216,p216);
        c216 = (a216^2-b216^2)^0.5;

        x216 = (a216)*cos(theta)-c216; % since a > b, focus (sun) is offset on the x axis
        y216 = (b216)*sin(theta);



        % graph

        theta = linspace(0,2*pi,500);
        sun_y = paray(rSun,theta); % sun y data
        sun_x = parax(rSun,theta);

        figure
        title(sprintf('Orbital Transfer From 16 Psyche to 216 Kleopatra'))
        hold on
        plot(x16,y16,'b','Linewidth',2) % 16 Psyche orbit blue
        plot(x216,y216,'g','Linewidth',2) % 216 Kleopatra orbit green
        plot(xTsfr,yTsfr,'r','Linewidth',2) % transfer orbit red
        plot(sun_x,sun_y,'y','Linewidth',2)
        
%         temp = linspace(1,rSun,length(1:1000000:rSun));
%         for i = 1:length(1:1000000:rSun) % plot sun as yellow
%             
%             sun_x = parax(temp(i),theta);
%             plot(sun_x,sun_y,'y')
%             
%         end
        
        
        xlabel('m')
        ylabel('m')
        set(gca, 'color', [0,0,0])
        axis equal
        lgd = legend('16 Psyche','216 Kleopatra','Transfer Orbit','Sun');
        %c = lgd.TextColor;
        lgd.TextColor = [1,1,1];
        hold off

    case 0 % fire same direction as asteroid
        
        % inner asteroid
        va16 = derivedV(muSun,a16,ra16);
        vp16 = derivedV(muSun,a16,rp16);
        p16 = semilatus1(a16,e16);
        
        
        % outer asteroid
        va216 = derivedV(muSun,a216,ra216);
        vp216 = derivedV(muSun,a216,rp216);
        p216 = semilatus1(a216,e216);

        % Hohman transfer orbit
        rpTsfr = rp16; % by nature of outer orbit geometry firing at apogee
        raTsfr = ra216; % by nature of outer orbit geometry firing at apogee
        
        eTsfr  = eccentricity1(raTsfr,rpTsfr);
        aTsfr = majorAxisHalf1(raTsfr,rpTsfr);
        pTsfr = semilatus1(aTsfr,eTsfr); 
        vaTsfr = derivedV(muSun,aTsfr,raTsfr);
        vpTsfr = derivedV(muSun,aTsfr,rpTsfr);

        tTsfr = TOF_ellipse(aTsfr,muSun)/2;
        t16 = TOF_ellipse(a16,muSun)/2;
      
        dv1 = abs(vp16-vpTsfr); % dv leaving psyche
        dv2 = abs(vaTsfr-va216); % dv to enter kleopatra


        theta = linspace(0,2*pi,500);

        bTsfr = bGeo(aTsfr,pTsfr);
        cTsfr = (aTsfr^2-bTsfr^2)^0.5;

        xTsfr = (aTsfr)*cos(theta)-cTsfr; % since a > b, focus (sun) is offset on the x axis
        yTsfr = (bTsfr)*sin(theta);

        b16 = bGeo(a16,p16);
        c16 = (a16^2-b16^2)^0.5;

        x16 = (a16)*cos(theta)-c16; % since a > b, focus (sun) is offset on the x axis
        y16 = (b16)*sin(theta);
        
        b216 = bGeo(a216,p216);
        c216 = (a216^2-b216^2)^0.5;

        x216 = (a216)*cos(theta)-c216; % since a > b, focus (sun) is offset on the x axis
        y216 = (b216)*sin(theta);



        % graph

        theta = linspace(0,2*pi,500);
        sun_y = paray(rSun,theta); % sun y data

        figure
        title(sprintf('Orbital Transfer From 16 Psyche to Outer Orbit'))
        hold on
        plot(x16,y16,'b','Linewidth',2) % 16 Psyche orbit blue
        plot(x216,y216,'g','Linewidth',2) % 216 Kleopatra orbit green
        plot(xTsfr,yTsfr,'r','Linewidth',2) % transfer orbit red

        temp = linspace(1,rSun,length(1:1000000:rSun));
        for i = 1:length(1:1000000:rSun) % plot sun as yellow
            
            sun_x = parax(temp(i),theta);
            plot(sun_x,sun_y,'y')
            
        end
       
        xlabel('m')
        ylabel('m')
        set(gca, 'color', [0,0,0])
        axis equal
        legend off
        hold off

end

end

function [dv4,tCap,rp,vp] = capture(mu216,muSun,a216,r216,vInf)
soi216 = SOI(mu216,muSun,a216);
b216 = impactParameter(vInf,mu216,r216);

% pick d! i chose half the distance between b and soi, but can be weighted
% accordingly
d = b216 + 0.50*(soi216-b216);

% find dv to enter circular orbit, graph resulting trajectories

% incoming flyby
E = energy2(vInf,mu216,soi216);
H = angularMomentum4(vInf,d);
p = semilatus2(H,mu216);
e  = eccentricity2(E,H,mu216);
a = majorAxisHalf3(mu216,E);
rp = radiusPeriapsis(a,e);
vp = velocity(E,mu216,rp);

% colinear firing
vcs = circularVelocity(mu216,rp);

dv4 = abs(vp-vcs);

% TOF

nu = trueAnomoly(e,p,soi216);
tCap = TOF_hyp(nu,mu216,e,a);

% graph

theta = linspace(0,2*pi,500);
[dx,dy] = para(rp,theta); % orbit D x and y
psy_y = paray(r216,theta); % psyche y data
[ex,ey] = paraHyperbola(a,p); % escape x and y

figure
title(sprintf('(Ideal) 216 Kleopatra Capture with Vinf = %g, circular park at %gm',vInf,rp))
hold on
plot(dx,dy,'b','Linewidth',2) % orbit D
for i = 1:100:r216 % plot kleopatra as cyan
    psy_x = parax(i,theta);
    plot(psy_x,psy_y,'c')
end
plot(ex,ey, 'r','Linewidth',2) % escape
xlabel('m')
ylabel('m')
set(gca, 'color', [0,0,0])
axis equal
hold off

end

function [dv4,tCap,rp,vp] = flyby(mu216,muSun,a216,r216,vInf)
soi216 = SOI(mu216,muSun,a216);
b216 = impactParameter(vInf,mu216,r216);

% pick d! i chose half the distance between b and soi, but can be weighted
% accordingly
d = b216 + 0.50*(soi216-b216);

% find dv to enter circular orbit, graph resulting trajectories

% incoming flyby
E = energy2(vInf,mu216,soi216);
H = angularMomentum4(vInf,d);
p = semilatus2(H,mu216);
e  = eccentricity2(E,H,mu216);
a = majorAxisHalf3(mu216,E);
rp = radiusPeriapsis(a,e);
vp = velocity(E,mu216,rp);

% colinear firing
vcs = circularVelocity(mu216,rp);

dv4 = abs(vp-vcs);

% TOF

nu = trueAnomoly(e,p,soi216);
tCap = TOF_hyp(nu,mu216,e,a);

% graph

theta = linspace(0,2*pi,500);
[dx,dy] = para(rp,theta); % orbit D x and y
psy_y = paray(r216,theta); % psyche y data
[ex,ey] = paraHyperbola(a,p); % escape x and y

figure
title(sprintf('216 Kleopatra Flyby with Vinf = %g m/s, closest pass at %gm',vInf,rp))
hold on
plot(dx,dy,'b','Linewidth',2) % orbit D
for i = 1:100:r216 % plot kleopatra as cyan
    psy_x = parax(i,theta);
    plot(psy_x,psy_y,'c')
end
plot(ex,ey, 'r','Linewidth',2) % escape
xlabel('m')
ylabel('m')
set(gca, 'color', [0,0,0])
axis equal
hold off

end