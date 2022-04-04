%% Clean up the workspace, and define constants, and things that we will care about for all versions
clear
close all
clc


%  Define the constants per the problem.  We use I_eq as 6 as this will
%  allow us to have a bistable system
C = 10;
gl = 19;
El = -67;
gna = 74;
Ena = 60;
V_half = 1.5;
k = 16;
I_eq = 6;

%  Define the timestep spacing. Euler's method has failed!
delta_t = 0.3;

%  Define the initial condition(s).
v_values = 0;

%  Define the timespan
tspan = 0:delta_t:20;
num_time_steps = length(tspan);

%  Use the boltzman equation (really, a sigmoid) to define h_infinity. Also
%  define its derivative
m_inf = @(V) 1./(1+exp((V_half - V)./k));
dm = @(V) -((-1./k).*(exp((V_half - V)./k)))./((1+exp((V_half - V)./k)).^2);
vdot = @(t,V) (I_eq - (gl.*(V-El) + gna.*m_inf(V).*(V-Ena)))/C;

%% This evolves the solution over time using Euler's method
Y1 = nan(num_time_steps, length(v_values)); 
Y1(1) = v_values;
T1 = nan(num_time_steps,1);
T1(1) = 0;
for euler_counter = 1:num_time_steps-1
    T1(euler_counter+1) = euler_counter*delta_t;
    Y1(euler_counter+1) = Y1(euler_counter) + vdot(T1(euler_counter),Y1(euler_counter))*delta_t;
end

figure(1)
plot(T1,Y1, '-o', 'linewidth', 3)
xlabel('time')
ylabel('V')
title('Euler''s method')
h = gca;
h.FontSize = 15;


%% This evolves the solution over time using RK4 with anonymous functions

Y2 = nan(num_time_steps, length(v_values)); 
Y2(1) = v_values;
T2 = nan(num_time_steps,1);
T2(1) = 0;
for rk_counter = 1:num_time_steps-1
    T2(rk_counter+1) = rk_counter*delta_t;
    Y2(rk_counter+1) = rk4(Y2(rk_counter), delta_t, vdot, T2(rk_counter));
end

figure(2)
plot(T2,Y2, '-o', 'linewidth', 3)
xlabel('time')
ylabel('V')
title('RK4')
h = gca;
h.FontSize = 15;


%% This evolves the solution over time using ode45 with external functions
[T3, Y3] = ode45(@(t,y) vdot(t,y), tspan, v_values);


figure(3)
plot(T3,Y3, '-o', 'linewidth', 3)
xlabel('time')
ylabel('V')
title('ODE45 with regular functions')
h = gca;
h.FontSize = 15;



function result = rk4(Un, k, f, t)
%  Function to compute the update for one time step of rk4
%  USAGE
%  Inputs
%  Un - Solution at the current timestep
%  k - timestep spacking
%  f - anonymous function of the function that we are integrating
%  t - current timestep
%
%  Output
%  result - value at the next timestep

y1 = Un;
y2 = Un + k.*((1/2).*f(t, y1));
y3 = Un + k.*((1/2).*f(t+0.5*k,y2));
y4 = Un + k.*((1)*f(t+0.5*k,y3));

result = Un + k.*((1/6).*(f(t,y1))+...
                  (1/3).*(f(t+0.5*k,y2))+...
                  (1/3).*(f(t+0.5*k,y3))+...
                  (1/6).*(f(t+k,y4)));
end
