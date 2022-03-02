function Example_10_6
% 
%
% This function solves Example 10.6 by using MATLAB’s ode45 to numerically
% integrate Equations 10.89 (the Gauss planetary equations) to determine
% the J2 perturbation of the orbital elements.
%
% User M-functions required: None
% User subfunctions required: rates
% ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
%...Preliminaries:
close all; clear all; clc
%...Conversion factors:
hours = 3600; %Hours to seconds
days = 24*hours; %Days to seconds
deg = pi/180; %Degrees to radians
%...Constants:
mu = 42646; %Gravitational parameter (km^3/s^2)
RE = 3389.5; %Earth's/Mars' radius (km)
J2 = 1955.5e-6; %Earth/s/Mars' J2
%...Initial orbital parameters (given):
rp0 = RE + 500; %perigee radius (km)
ra0 = RE + 550; %apogee radius (km)
RA0 = 45*deg; %Right ascencion of the node (radians)
i0 = 63*deg; %Inclination (radians)
w0 = 270*deg; %Argument of perigee (radians)
TA0 = 40*deg; %True anomaly (radians)
%...Initial orbital parameters (inferred):
e0 = (ra0 - rp0)/(ra0 + rp0); %eccentricity
h0 = sqrt(rp0*mu*(1 + e0)); %angular momentum (km^2/s)
a0 = (rp0 + ra0)/2; %Semimajor axis (km)
T0 = 2*pi/sqrt(mu)*a0^1.5; %Period (s)
%...Store initial orbital elements (from above) in the vector coe0:
coe0 = [h0 e0 RA0 i0 w0 TA0];
%...Use ODE45 to integrate the Gauss variational equations (Equations
% 12.89) from t0 to tf:
t0 = 0;
tf = 700*days;
nout = 5000; %Number of solution points to output for plotting purposes, inizialmente 5000
tspan = linspace(t0, tf, nout);
options = odeset(...
'reltol', 1.e-8, ...
'abstol', 1.e-8, ...
'initialstep', T0/1000);

sum = 0;
sum1 = 0;
[rsat0, v0] = sv_from_coe(coe0,mu);
r0 = norm(rsat0)
tx = [(3/2)*(J2*mu*(RE^2)/(r0^4))*(rsat0(1)/r0)*(5*(rsat0(3)^2)/(r0^2)-1) (3/2)*(J2*mu*(RE^2)/(r0^4))*(rsat0(2)/r0)*(5*(rsat0(3)^2)/(r0^2)-1) (3/2)*(J2*mu*(RE^2)/(r0^4))*(rsat0(3)/r0)*(5*(rsat0(3)^2)/(r0^2)-3)];
temp0 = norm(tx)
y0 = [coe0 temp0]';


[t,y] = ode45(@rates, tspan, y0, options);
%...Assign the time histories mnemonic variable names:
h = y(:,1);
e = y(:,2);
RA = y(:,3);
i = y(:,4);
w = y(:,5);
TA = y(:,6);
temp = y(:,7);
%...Plot the time histories of the osculatinig elements:
figure(1)
subplot(5,1,1)
plot(t/86400,(RA - RA0)/deg)
title('Right Ascension (degrees)')
xlabel('days')
grid on
grid minor
axis tight
subplot(5,1,2)
plot(t/86400,(w - w0)/deg)
title('Argument of Perigee (degrees)')
xlabel('days')
grid on
grid minor
axis tight
subplot(5,1,3)
plot(t/86400,h - h0)
title('Angular Momentum (km^2/s)')
xlabel('days')
grid on
grid minor
axis tight
subplot(5,1,4)
plot(t/86400,e - e0)
title('Eccentricity')
xlabel('days')
grid on
grid minor
axis tight
subplot(5,1,5)
plot(t/86400,(i - i0)/deg)
title('Inclination (degrees)')
xlabel('days')
grid on
grid minor
axis tight
%Plot the a_J2 effect
figure(2)
plot(t/86400, temp)
title('norm A_J2')
xlabel('days')
grid on
grid minor
axis tight
%...Subfunction:
% 
function dfdt = rates(t,f)
% 
%
% This function calculates the time rates of the orbital elements
% from Gauss’s variational equations (Equations 12.89).
% –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
%...The orbital elements at time t:
h = f(1);
e = f(2);
RA = f(3);
i = f(4);
w = f(5);
TA = f(6);
coe = [h e RA i w TA];
[rsat, v] = sv_from_coe(coe,mu);
r =h^2/mu/(1 + e*cos(TA)); %The radius
u = w + TA; %Argument of latitude
%...Orbital element rates at time t (Equations 12.89):
hdot = -3/2*J2*mu*RE^2/r^3*sin(i)^2*sin(2*u);
edot = ...
3/2*J2*mu*RE^2/h/r^3*(h^2/mu/r ...
*(sin(u)*sin(i)^2*(3*sin(TA)*sin(u) - 2*cos(TA)*cos(u)) - sin(TA)) ...
-sin(i)^2*sin(2*u)*(e + cos(TA)));
edot = 3/2*J2*mu*RE^2/h/r^3 ...
*(h^2/mu/r*sin(TA)*(3*sin(i)^2*sin(u)^2 - 1) ...
-sin(2*u)*sin(i)^2*((2+e*cos(TA))*cos(TA)+e));
TAdot = h/r^2 + 3/2*J2*mu*RE^2/e/h/r^3 ...
*(h^2/mu/r*cos(TA)*(3*sin(i)^2*sin(u)^2 - 1) ...
+ sin(2*u)*sin(i)^2*sin(TA)*(h^2/mu/r + 1));
RAdot = -3*J2*mu*RE^2/h/r^3*sin(u)^2*cos(i);
idot = -3/4*J2*mu*RE^2/h/r^3*sin(2*u)*sin(2*i);
wdot = 3/2*J2*mu*RE^2/e/h/r^3 ...
*(-h^2/mu/r*cos(TA)*(3*sin(i)^2*sin(u)^2 - 1) ...
- sin(2*u)*sin(i)^2*sin(TA)*(2 + e*cos(TA)) ...
+ 2*e*cos(i)^2*sin(u)^2);
a_grav = -(mu*rsat)/(r^3);
%a_tot = a_grav
a_J2 = [(3/2)*(J2*mu*(RE^2)/(r^4))*(rsat(1)/r)*(5*(rsat(3)^2)/(r^2)-1) (3/2)*(J2*mu*(RE^2)/(r^4))*(rsat(2)/r)*(5*(rsat(3)^2)/(r^2)-1) (3/2)*(J2*mu*(RE^2)/(r^4))*(rsat(3)/r)*(5*(rsat(3)^2)/(r^2)-3)];
norm(a_J2);  %km/s^2
temp = norm(a_J2);
sum = sum + norm(a_J2);
sum1 = sum1 + 1;
%...Pass these rates back to ODE45 in the array dfdt:
dfdt = [hdot edot RAdot idot wdot TAdot temp]';
end %rates

a_J2_average = sum/sum1 
Delta_v_J2 = a_J2_average*(60*60*24)  % km/s

% 
end %Example_10_6

