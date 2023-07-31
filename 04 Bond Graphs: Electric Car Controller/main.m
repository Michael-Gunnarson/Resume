% Group 10
% Lab 5
% EME 171 Fall 2022

close all
clear
clc

%% Variables

global Rw Lw Tm M bt R Gr Cr g Cd rho Af T1 Kp Ki swap dt

% system constants:
Rw = 0.3; % ohms, armature winding resistance
Lw = 0.015; % henry, armature winding inductance
Tm = 1.718; % Wber, transduction coefficient (what???)
M = 2200; % kg, vehicle mass
bt = 0.05; % N*m*s/rad, drive shaft friction
R = 0.2; % M, wheel radius
Gr = 5/1; % Gear ratio
Cr = 0.006; % rolling resistance coefficient
g = 9.81; % m/s^2, gravity
Cd = 0.32; % drag coefficient
rho = 1.21; % kg/m^3, air density
Af = 2.05; % m^2, vehicle frontal area


%% Initial Conditions

% we have two states:
% p3: from inductance, p9: from car mass
% since we wnat to see the distance travelled by the car, we will
% also integrate f9

p3_0 = 0;
p9_0 = 0;
d_0 = 0; % not a state but integration is desired
dref_0 = 0;

power_in_0 = 0;

initials = [p3_0,p9_0,d_0,dref_0,power_in_0];


%% Part 1: Modelling and Implementation: Step Input Test

% timesteps:
T1 = 5;
dt = 0.01;
tspan = 0:dt:(T1+4); % simulate 4 seconds after step

swap = "Part 1";
[t,s] = ode45(@eqns5,tspan,initials);

figure
plot(t,s(:,2)/M), grid on
title('Car Velocity for Step Input of 100 Volts at T=5 sec')
ylabel('velocity (m/s)')
xlabel('time (s)')


%% Part 2: PI Controller Design

% intuition: Kp affects speed of response, Ki affects steady state error
% from experimentation: Kp ^ means rise time goes down, but settling time
% and overshoot go up

Ki = 15;
Kp = 15;

% initialize:
overshoot = 1;
settle = 5;
rise = 1;

% while overshoot  .10 || (overshoot == 'N/A')
T1 = 1;
Tend = 25;
dt = 0.01;
tspan = 0:dt:Tend;
while settle > 2.0 || rise > 0.5 || any(overshoot > .10)
   if all(overshoot == 'N/A')
       Kp = Kp + 1;
       Ki = Ki + 1;
   swap = "Part 2";
   [t,s] = ode45(@eqns5,tspan,initials);
   [rise,settle,overshoot] = find_response(t,s);
   end

while settle > 2.0
  
   Ki = Ki + 1;

   swap = "Part 2";

   [t,s] = ode45(@eqns5,tspan,initials);
   [rise,settle,overshoot] = find_response(t,s);
end

while rise > 0.5
   Kp = Kp + 1;

   swap = "Part 2";
   [t,s] = ode45(@eqns5,tspan,initials);
   [rise,settle,overshoot] = find_response(t,s);
end

while any(overshoot > .10)
   if all(overshoot == 'N/A')
       break
   end
   Kp = Kp -1;
   Ki = Ki -1;

   swap = "Part 2";
   [t,s] = ode45(@eqns5,tspan,initials);
   [rise,settle,overshoot] = find_response(t,s);
 
end

if settle < 2.0 || rise < 0.5 || any(overshoot < .10)
tspan = 0:dt:(T1+4);
% disp('we made it to this part of the loop')

% I could not tell you why, but the exact breakdown of tspan has huge
% implications on the response of the controller.  This is such a janky
% solution, I just switched it to the tspan I want to use when we start
% getting viable solutions and solve it again.  This is the only way to get
% my tuning response to match my experimental responses near the critical
% result

end
end


%% print controller response characteristics

tspan = 0:dt:(T1+4); % four seconds after vref starts
swap = "Part 2";
[t,s] = ode45(@eqns5,tspan,initials);
[rise,settle,overshoot] = find_response(t,s)
sprintf("Ki = %g",Ki)
sprintf("Kp = %g",Kp)


%% Plot Vref

figure
plot(t,s(:,2)/M), grid on
title('Velocity Response for Vref = 1 at T = 1 sec (step input)')
ylabel('velocity (m/s)')
xlabel('time (s)')


%% Part 3: Testing

dt = 0.01;
tspan = 0:dt:300;
swap = "Part 3";
[t,s] = ode45(@eqns5,tspan,initials);

vref = zeros(length(t),1);
for i = 1:length(t)
   vref(i) = LA92Oracle(t(i));
   [ds(i,:) ext(i,:)] = eqns5(t(i), s(i,:));
end

uin = ext(:,1);
p9_dot = ext(:,2);
iin = ext(:,3);
vout = ext(:,4);
d = s(:,4);
power_integral = s(:,5);

Pin = uin.*iin;
Pout = p9_dot.*vout;

Pin_acc = Pin(Pout>0);
Pout_acc = Pout(Pout>0);
efficiency = mean(Pout_acc)/mean(Pin_acc)
meters_over_joules = max(d)/max(power_integral)

figure
plot(t,s(:,2)/M), grid on
title('Autonomous vehicle velocity response for LA92 velocity profile')
ylabel('velocity (m/s)')
xlabel('time (s)')

figure
hold on
plot(t,s(:,2)/M), grid on
plot(t,vref)
title('Autonomous vehicle velocity compared to LA92 velocity profile')
legend('Autonomous Electric Vehicle','LA92 reference')
ylabel('velocity (m/s)')
xlabel('time (s)')
axis([32,54,4,8])
hold off


function [rise,settle,overshoot] = find_response(t,s)
global Rw Lw Tm M bt R Gr Cr g Cd rho Af T1 Kp Ki dt swap
% analyze response characteristics: easily solvable since we know vref
   v = s(:,2)/M;
   % find indexes:

   % rise time is time between 10% and 90% of final value. vref = 1
   rise_i = find(v >= .1,1,"first");
   rise_f = find(v >= .9,1,"first");

   rise = t(rise_f) - t(rise_i);

   % to find settling: if dy/dx < 0, find where, after max, v crosses 1 + .02
   % [figure out how to bound equation after this point? aka it has too many oscillations]
   % if dy/dx > 0 for entire system, find where v crosses 1 - .02

   % numerically approximate derivative at every midpoint
   v2 = v(2:end);
   v1 = v(1:end-1);
   t2 = t(2:end);
   t1 = t(1:end-1);

   dvdt = (v2-v1)./(t2-t1);


   if any(dvdt < 0) % if dv/dt < 0, response overshoots
       over = find(v >= (1 +.02)); % find all above error response
       under = find(v <= (1 - .02)); % find all below error response
   
       set_index = max([over;under]);
       %settle = t(set_index + 1) - T1; % add one to index to get into the acceptable error region
       settle = t(set_index) - T1;
   
   else % if no dv/dt < 0, system is overdamped and does not cross infinite value
       set_index = find(v >= (1 -.02), 1,"first"); % approach from bottom
       settle = t(set_index);
   end
 
   
   
   if any(v > 1) % check if velocity response overshoots
       overshoot = max(v) - 1;
   else
       overshoot = 'N/A';
   end
end
