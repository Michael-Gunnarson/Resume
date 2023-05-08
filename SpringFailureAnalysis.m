% EME 150A Project 2 
% Team Morpheus
% Code written by Michael Gunnarson 
% with the help of Adela Alcazar

close all
clear
clc

%% Material Properties

% matweb.com for
% 17-7 PH Stainless Steel, CH900, wire
Sut = abs(229900+364800)/2; % psi, ultimate tensile strength

% Pulled from 17-7 PH Stainless Steel CH900 plate sheet or strip
Sy = 231000; % psi, Yield strength (needed for goodman line)
E = 29600e3; % psi 


% empirically
if(Sut > 200e3) % psi
    Se = 100e3; % psi
elseif(Sut <= 200e3)
    Se = 0.5*Sut; % psi
end




%% Fixed Parameters

% if you want to change geometries, change C and safety factor
theta = deg2rad(120); % angle for one full cycle
C = 12; % Spring index


%% Iterate Geometry

d = 0.14; % inches, wire diameter initial guess
N = 4; % number of active coils, to be iterated later
sf = 0.1; % safety factor initialization


% loop d to get desired safety factor
iter = 1;
while sf < 1.5 && iter < 100000

    % d iterated, C constant ==> changes D.  Find new T
    D = find_D(C,d); % mean coil diameter, inches
    T = find_desired_T(D); % since we want 60 lbf to be the highest force 
                           % at 1.5 inches from the top of the bar, 
                           % we use this function and the radius of mean diameter
                           % to find the desired torque in lbf*in

    % sigma and safety factor
    Sigma_max = sigma_max(T,d,C); % psi
    Sigma_min = 0; % designed with no reversability in mind
    Sm = 0.5*(Sigma_max+Sigma_min); % mean stress
    Sa = 0.5*(Sigma_max-Sigma_min); % amplitude stress
    
    sf = SF(Sut,Se,Sm,Sa); % safety factor

    % change value till convergence 
    d = d + 0.01;
    iter = iter + 1;
end

d = d - 0.01; % reset to value that actually broke loop

if iter == 100000
    disp("you're wrong")
end



Fmax = 100; % get the loop started
% loop N to get desired force distribution (forces desired k)
iter = 1;
while Fmax > 60 &&  iter < 100000
    
    % find 
    k = spring_const(E,d,C,N); % lbf*in
    f = force_distribution(k,theta,D); % lbf at each location along bar
                                       % assuming average hand width 3
                                       % inches we have 4 positions along
                                       % beam without overlap to be marked
                                       % by manufacturer
    Fmax = max(f);
    %disp(f)
    
    N = N + 0.1;
    iter = iter + 1;
end

N = N - 0.1; % reset to value that broke loop
N = N - 0.1; % reset to one iteration before, to go just over 60 lbf max force

% re-calculate with the correct N value
k = spring_const(E,d,C,N);
f = force_distribution(k,theta,D);
Fmax = max(f);

if iter == 100000
    disp("you're wrong")
end



%% Display Wire Geometry Parameters

disp('Geometric Parameters:')
sprintf('Mean Coil Diameter = %g in',D)
sprintf('Wire Diameter = %g in',d)
sprintf('Spring Index = %g',C)
sprintf('Number of Coils = %g',N)
sprintf('Maximum Stress = %g psi',Sigma_max)
sprintf('Safety Factor = %g',sf)
sprintf('At the four locatoins on the bar, the forces needed to fully bend the device are as follows (assuming average hand width about 3 inches to find resultant force placements:)')
L = [1.5,4.5,7.5,10.5]; % handle distribution
for i = 1:4
    sprintf('At %g inches from the top of the bar, you need %g lbf to move the handles 120 degrees',L(i),f(i))
end


%% Print Fatigue Diagram

title = "Fatigue Diagram, as designed";
fatigue_diagram(Se,Sy,Sut,Sa,Sm,title)
sprintf('Location on fatigue diagram shows us we have infinite fatigue life!')


%% One Last Calculation: Bending the wrong direction

% leave spring constant k as designed to see how it fairs if used
% incorrectly

Sigma_max_ = sigma_max(T,d,C); % psi, bent in correct direction
% if the design were to bend in the opposite direction until horizontal: 
theta_ = deg2rad(120-180);
T_ = find_T(k,theta_); % used not as designed
Sigma_min_ = sigma_max(T_,d,C);
Sm_ = 0.5*(Sigma_max_+Sigma_min_); % mean stress
Sa_ = 0.5*(Sigma_max_-Sigma_min_); % amplitude stress
    
sf_horz = SF(Sut,Se,Sm_,Sa_); % safety factor
disp('If used incorrectly, pulling bars to horizontal configuration:')
sprintf('Safety Factor for horizontal configuration (used incorrectly) = %g',sf_horz)
sprintf('Minimum Stress = %g psi',Sigma_min_)

title = 'Fatigue Diagram used improperly, horizontal configuration';
fatigue_diagram(Se,Sy,Sut,Sa_,Sm_,title)



% if the design were to bend to zero degrees, passing the horizontal
% component on the way:

theta_ = deg2rad((120-180) - 180);

T_ = find_T(k,theta_); % used not as designed
Sigma_min_ = sigma_max(T_,d,C);
Sm_ = 0.5*(Sigma_max_+Sigma_min_); % mean stress
Sa_ = 0.5*(Sigma_max_-Sigma_min_); % amplitude stress
    
sf_zero = SF(Sut,Se,Sm_,Sa_); % safety factor
disp('If used incorrectly, pulling bars past horizontal configuration all the way to parallel:')
sprintf('Safety Factor for parallel bars in wrong direction (used incorrectly) = %g',sf_zero)
sprintf('Minimum Stress = %g psi',Sigma_min_)

title = 'Fatigue Diagram used improperly, parallel bars in wrong direction';
fatigue_diagram(Se,Sy,Sut,Sa_,Sm_,title)


%% Functions

function D = find_D(C,d)
D = C*d;
end

function T = find_desired_T(D)
% T = 2FL
% at 1.5 inches down handle, desired force is 60 lbf on each side

T = 2*60*(1.5+D/2);
end

function s = sigma_max(T,d,C)
K = big_K(C);
s = K*32*T/(pi*d^3);
end

function k = big_K(C)
k = (4*C^2-C-1)/(4*C*(C-1));
end

function n = SF(Sut,Se,Sm,Sa)
n = (Sut*Se)/(Sm*Se+Sa*Sut);
end

function k = spring_const(E,d,C,N)
k = E*d^3/(64*C*N);
end

function f = force_distribution(k,theta,D)

L = [1.5,4.5,7.5,10.5]; % handle distribution
L = L + D/2; % including distance from spring

% for torsional spring, torque response is constant
% T = k*theta
% T = 2*F*L (for our specific machine, force applied on both sides),
% then F = k*theta/(2*L)

f = k*theta./(2*L);

end

function fatigue_diagram(Se,Sy,Sut,Sa,Sm,Title)

% convert from psi to kpsi
a = 1e3;
Se = Se/a;
Sy = Sy/a;
Sut = Sut/a;
Sa = Sa/a;
Sm = Sm/a;

% static yielding lines
x1 = -Sy:0;
x2 = 0:Sy;

m1 = (Sy-0)/(0-(-Sy));
m2 = (0-Sy)/(Sy-0);

b1 = Sy;
b2 = Sy;

y1 = m1*x1 + b1;
y2 = m2*x2 + b2;


% goodman line
% point of intersection/horizontal component:
x_int = (Se-b1)/m1;
x_horz = x_int:0;
y_horz = Se*ones(1,length(x_horz));

m_good = (0-Se)/(Sut-0);
b_good = Se;
x_good = 0:Sut;
y_good = m_good*x_good+b_good;

figure
title(Title)
hold on
plot(x1,y1, 'black')
plot(x2,y2, 'black')
plot(x_horz,y_horz, 'red')
plot(x_good,y_good, 'red')

scatter(Sm,Sa,'x','blue')
xlabel('Mean Stress (kpsi)')
ylabel('Stress Amplitude (kpsi)')
legend('Static Failure Line','Static Failure Line','Goodman Horizontal Component','Goodman Line','Part Location')
hold off
end

function T = find_T(k,theta)
% theta must be in radians
T = k*theta;
end