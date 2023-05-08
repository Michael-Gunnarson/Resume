% Team Morpheus
% code written by Michael Gunnarson
% 10/23/2022
% Written in MATLAB R2022a

clear
clc

% This matlab script analyzes various belt forces within an amazon record
% player.  The initial force is calculated using geometry, a centrifugal
% force calculated using mass of the belt and speed of rotation, and the
% torque from the spec sheet is added to either side of the equation in
% order to find the tight side and loose side tensions.  The coefficient of
% friction is then found numerically.

%% Initial Measurements

D = 5.1; % in, plate diameter
d = 4.5; % mm, motor diameter
c = 3.1; % in, distance between axles

L_orig = 8; % in, doubled up on itself
L_orig = L_orig*2;

% convert to meters
D = in_m(D);
d = mm_m(d);
c = in_m(c);
L_orig = in_m(L_orig);


%% Test to find E
mass = 187; % g, pocket knife
weight = mass/1000 * 9.81;
% since we've doubled it up, the force felt by each "member" is F/2
weight = weight/2;

% rubber band measurements
b = 1.5; % mm
h = 0.5; % cm

b = mm_m(b);
h = cm_m(h); % convert to meters

band_area = b*h;

% after stretching
L_test = 8.6; % in, doubled up on itself
L_test = L_test*2; % unwinding the band
L_test = in_m(L_test);

delta_test = (L_test - L_orig);

% epsilon_test = delta_test/L_orig

E = weight*L_orig/(band_area*delta_test); % Pa, end of test
fprintf('The Experimental Elastic Modulus is %g Pa \n',E)

%% Calculate Forces Within Band
L_stretch = length_stretch(D,d,c);
delta = (L_stretch-L_orig);

epsilon = delta/L_orig;

Fi = band_area*E*epsilon;
fprintf('The initial tension in stationary band is %g N \n',Fi)


% solve for fc (belting equation)
% mass of belt
mass = 2/1000; % kg
r = d/2;

rpm = [33,45,78];
omega_follow = rpm*(2*pi)/60; % rad/s

omega_drive = D/d*omega_follow;
Fc = mass*(r*omega_drive).^2;

for i = 1:3
fprintf('The Centrifugal force due to moving belt is %g N ',Fc(i))
fprintf('for %i rpm \n',rpm(i))
end

% from Spec sheet
% https://www.consumer-parts.com/product.sc?productId=29486
T = 8; % gram force centimeters, torque

% should be 0.000784532
T = T*9.80665e-5;


F1 = zeros(3,1);
F2 = F1;
for i = 1:3
F1(i) = Fi + Fc(i) + T/d;
F2(i) = Fi + Fc(i) - T/d;
end


for i = 1:3
    fprintf('The tight side tension is %g and the loose side tension is %g N ',F1(i),F2(i))
    fprintf('for %i rpm \n',rpm(i))
end
fprintf('Tight side and loose side tension arises from the rotating motor and coefficient of friction')


% find where g is equal to zero
fid = fi_d(D,d,c);
g = @(f) T/d*(exp(f*fid)+1)/(exp(f*fid)-1) - Fi;

f = Secant_Method(g,0.1,0.01);

fprintf('The Coefficient of Friction found Numerically is %g', f)


%% Functions

% translate to meters:
function m = mm_m(mm)
m = mm/1000;
end

function m = in_m(in)
m = in*2.54/100;
end

function in = m_in(m)
in = m*100/2.54;
end

function m = cm_m(cm)
m = cm/100;
end


% change in length based on geometry
function fiD = fi_D(D,d,c)
fiD = pi + 2*asin((D-d)/(2*c));
end

function fid = fi_d(D,d,c)
fid = pi - 2*asin((D-d)/(2*c));
end

function L_prime = length_stretch(D,d,c)
    fiD = fi_D(D,d,c);
    fid = fi_d(D,d,c);
    L_prime = (4*c^2-(D-d)^2)^0.5 + 0.5*(D*fiD + d*fid);
end


% initial force, Newton's method, etc

function xNew = Secant_Method(f,xn,xnminus1)
% the secan method is an offshoot of Newton's method, using the finite
% difference as an approximation for the derivative.  Requires two initial
% guesses in order to begin working

tol = 10^-5;
res = 1;
iter = 1;
%xi = xn;
%xiplus = xnminus1;

    %tStart = tic;
    while res > tol && iter < 10000
        
        df = (f(xn)-f(xnminus1))/(xn-xnminus1); 
        
xNew = xn - f(xn)/(df);

res = abs(f(xNew));  % because resolution for linear systems
                                     % is Ax-b, which should be zero, the
                                     % translation to the function case is
                                     % finding the difference between the
                                     % function out put and zero (y-0), or just
                                     % looking at the function output
        cor = correction(xNew,xn);
        
        iter = iter + 1;
        resolution_plot(iter) = res;
        correction_plot(iter) = cor;
        
        % shift everything to the right
        xnminus1 = xn;
        xn = xNew;
        
        
    end
    
%     tEnd = toc(tStart);
%    
%     figure
%     loglog(1:iter,resolution_plot)
%     xlabel('iteration')
%     ylabel('resolution')
%     title(sprintf('Secant method for initial guesses x=%G, x=%G',xi,xiplus))
%     
%     figure
%     loglog(1:iter,correction_plot)
%     xlabel('iteration')
%     ylabel('correction')
%     title(sprintf('Secant method for initial guesses x=%G, x=%G',xi,xiplus))
%  
%     tEnd
%     iter
%     res
    
end


function cor = correction(x1,x2)
    % x1 new
    % x2 old
    
   cor = max(max(abs(x1-x2)));

end
